"""Readers of RH's binary outputs (standard library only). The formats live
here. Arrays of doubles come back as array.array('d') in native byte order
(numpy: np.frombuffer(a) without a copy; they index, slice, iterate and
compare like tuples)."""
import struct
import sys
from array import array


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def doubles(raw, offset, n, big_endian=True):
    """n doubles of the bytes raw at offset (big-endian XDR, or native), as an
    array('d') in native byte order. struct.error, with the message of
    struct.unpack, if raw is too short."""
    chunk = raw[offset:offset + 8 * n]
    if len(chunk) != 8 * n:
        raise struct.error(f"unpack requires a buffer of {8 * n} bytes")
    a = array("d")
    a.frombytes(chunk)
    if big_endian and sys.byteorder == "little":
        a.byteswap()
    return a


def read_doubles(path):
    """A file of native-endian doubles (J.dat, Jgas.dat)."""
    raw = read_bytes(path)
    n = len(raw) // 8
    if len(raw) != 8 * n:
        raise struct.error(f"unpack requires a buffer of {8 * n} bytes")
    return doubles(raw, 0, n, big_endian=False)


def read_pops_array(path, lte=False):
    """(Nlevel, Nspace, NLTE populations as a flat list [level][depth]) of a
    pops.*.out file (the layout analyze_run.read_pops detects): XDR, an
    unknown header, then Nlevel and Nspace as ints, the populations and the
    LTE populations as doubles. With lte, the third item is the LTE
    populations (written after the NLTE ones)."""
    b = read_bytes(path)
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            start = off + 8 + (8 * nl * ns if lte else 0)
            return nl, ns, doubles(b, start, nl * ns)
    raise ValueError(f"{path}: cannot parse as a pops file (no layout with Nlevel, Nspace "
                     f"and 2 x Nlevel x Nspace doubles fits the size {len(b)})")


def read_pops(path):
    """NLTE populations of a pops.*.out file as a flat list."""
    return read_pops_array(path)[2]


def xdr_counted_string(b, off):
    """(string, next offset) of an XDR counted string at off in the bytes b,
    as RH writes them: xdr_short with the length (0: empty string, nothing
    more), then xdr_string (4-byte length, the bytes padded to 4)."""
    n = struct.unpack(">i", b[off:off + 4])[0]
    if n == 0:
        return "", off + 4
    m = struct.unpack(">i", b[off + 4:off + 8])[0]
    s = b[off + 8:off + 8 + m].decode()
    return s, off + 8 + 4 * ((m + 3) // 4)


def read_atmos_out(path, nspace):
    """(T, nHtot, abundances) of atmos.out (writeAtmos: ints Nlevel, Nelem,
    bool moving; T, ne, vturb; H populations [level][depth]; counted
    atmosphere ID; Nelem x (counted ID, weight, abundance)). nHtot is the sum
    of the H populations per depth point; abundances is a dict by upper-case
    element name."""
    b = read_bytes(path)
    nlev, nelem, _ = struct.unpack(">3i", b[:12])
    T = struct.unpack(f">{nspace}d", b[12:12 + 8 * nspace])
    off = 12 + 3 * 8 * nspace
    h = struct.unpack(f">{nlev * nspace}d", b[off:off + 8 * nlev * nspace])
    off += 8 * nlev * nspace
    _, off = xdr_counted_string(b, off)
    ab = {}
    for _ in range(nelem):
        el, off = xdr_counted_string(b, off)
        w, abund = struct.unpack(">2d", b[off:off + 16])
        off += 16
        ab[el.strip().upper()] = abund
    nh = [sum(h[l * nspace + k] for l in range(nlev)) for k in range(nspace)]
    return T, nh, ab


def read_first_molecule(path):
    """(ID, Tmin, Tmax) of the first molecule of molecules.out. Header:
    Nmolecule; counted ID (xdr short, xdr string); Nv, NJ; Ediss, Tmin,
    Tmax. ValueError if the header does not parse."""
    b = read_bytes(path)
    try:
        nmol = struct.unpack(">i", b[:4])[0]
        if nmol <= 0:
            raise ValueError(f"Nmolecule {nmol}")
        ident, off = xdr_counted_string(b, 4)
        if not ident:
            raise ValueError("empty molecule ID")
        ident = ident.strip()
        tmin, tmax = struct.unpack(">2d", b[off + 16:off + 32])    # after Nv, NJ, Ediss
    except (struct.error, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"{path}: no molecule header ({e})") from None
    return ident, tmin, tmax


def read_spectrum(path, nrays):
    """spectrum.out of an unpolarized run: Nspect (int), lambda, I
    [Nspect][nrays], big-endian XDR. Returns (lambda, rows of nrays values)."""
    raw = read_bytes(path)
    ns = struct.unpack(">i", raw[:4])[0]
    lam = doubles(raw, 4, ns)
    pos = 4 + 8 * ns
    I = doubles(raw, pos, ns * nrays)
    return lam, [I[i * nrays:(i + 1) * nrays] for i in range(ns)]


def read_stokes_spectrum(path, nrays):
    """spectrum.out of a polarized run: Nspect (int), lambda, I [Nspect][nrays]
    (doubles), a 12-byte record boundary, Q, U, V. All big-endian XDR.
    Returns (lambda, [I, Q, U, V]) with each Stokes array as a list of rows
    [Nspect][nrays]."""
    raw = read_bytes(path)
    ns = struct.unpack(">i", raw[:4])[0]
    lam = doubles(raw, 4, ns)
    size, pos = ns * nrays, 4 + 8 * ns
    I = doubles(raw, pos, size)
    pos += 8 * size + 12
    QUV = [doubles(raw, pos + 8 * size * i, size)
           for i in range(3)]
    if pos + 24 * size != len(raw):
        raise ValueError(f"{path}: unexpected size {len(raw)} for Nspect = {ns}, "
                         f"{nrays} rays (expected {pos + 24 * size})")
    return lam, [[q[l * nrays:(l + 1) * nrays] for l in range(ns)] for q in [I] + QUV]


def read_solveray_intensity(path, m):
    """The m intensities of an unpolarized solveray output (spectrum_*): mux,
    mu (2 doubles), I, the int Nspect = 0 of ray.input, big-endian XDR."""
    raw = read_bytes(path)
    if len(raw) != 16 + 8 * m + 4:
        raise ValueError(f"{path}: size {len(raw)}, expected {16 + 8 * m + 4} for {m} values")
    return doubles(raw, 16, m)


def read_solveray(path):
    """I, Q, U, V of a solveray output (spectrum_*): mux, mu (2 doubles), I,
    the int Nspect = 0 of ray.input, then Q, U, V, all big-endian XDR."""
    b = read_bytes(path)
    M = (len(b) - 20) // 32
    if 32 * M + 20 != len(b):
        raise ValueError(f"{path}: size {len(b)} does not fit the solveray output layout "
                         f"(32 M + 20 bytes)")
    off = 16 + 8 * M + 4
    return [doubles(b, o, M)
            for o in (16, off, off + 8 * M, off + 16 * M)]


VACUUM_TO_AIR_LIMIT = 200.0     # spectrum.h


def read_spectrum_file(path, stokes=None):
    """spectrum.out of either layout (writespect_xdr.c, big-endian XDR: int
    Nspect, double lambda[Nspect], double I[Nspect * Nrays * Nspace], a bool
    and VACUUM_TO_AIR_LIMIT, then Q, U and V of the size of I if the run was
    polarized). The size does not tell the layouts apart (Q, U, V read as
    part of I still make a multiple of Nspect values, Q-21): the record that
    must follow I (a bool, 0 or 1, and 200.0) does. stokes=True or False
    demands that layout, None finds it. Returns (raw bytes, lambda, arrays,
    Nspect) with arrays a dict {"I": [...]} or {"I", "Q", "U", "V"} of flat
    array('d') with the wavelength as the slowest index. ValueError (or
    struct.error for a file shorter than 4 bytes) if the file fits no
    layout."""
    raw = read_bytes(path)
    ns = struct.unpack(">i", raw[:4])[0]
    if ns <= 0 or len(raw) < 4 + 8 * ns + 12:
        raise ValueError(f"{path}: not a spectrum.out (Nspect={ns}, "
                         f"size {len(raw)})")
    lam = doubles(raw, 4, ns)
    body = len(raw) - 4 - 8 * ns - 12
    fits = []
    for narr in ((1, 4) if stokes is None else ((4,) if stokes else (1,))):
        if body <= 0 or body % (8 * narr) or (body // (8 * narr)) % ns:
            continue
        n_i = body // (8 * narr)
        marker = 4 + 8 * ns + 8 * n_i
        flag, limit = struct.unpack(">id", raw[marker:marker + 12])
        if flag in (0, 1) and limit == VACUUM_TO_AIR_LIMIT:
            fits.append((narr, n_i))
    if len(fits) != 1:
        want = {None: "either layout", True: "the polarized layout",
                False: "the unpolarized layout"}[stokes]
        raise ValueError(f"{path}: size {len(raw)} and contents do not fit "
                         f"{want} of spectrum.out (Nspect={ns})")
    narr, n_i = fits[0]
    off = 4 + 8 * ns
    arrays = {}
    for k, name in enumerate("IQUV"[:narr]):
        start = off if k == 0 else off + 8 * n_i + 12 + 8 * n_i * (k - 1)
        arrays[name] = doubles(raw, start, n_i)
    return raw, lam, arrays, ns
