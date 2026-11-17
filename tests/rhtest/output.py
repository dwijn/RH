"""Readers of RH's binary outputs (struct only). The formats live here."""
import struct


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def read_doubles(path):
    """A file of native-endian doubles (J.dat, Jgas.dat)."""
    raw = read_bytes(path)
    return struct.unpack(f"={len(raw) // 8}d", raw)


def read_pops_array(path):
    """(Nlevel, Nspace, NLTE populations as a flat list [level][depth]) of a
    pops.*.out file (the layout analyze_run.read_pops detects): XDR, an
    unknown header, then Nlevel and Nspace as ints, the populations and the
    LTE populations as doubles."""
    b = read_bytes(path)
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return nl, ns, struct.unpack(f">{nl * ns}d", b[off + 8:off + 8 + 8 * nl * ns])
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


def read_molecules(path, nspace):
    """(molecules, nHmin) of molecules.out (xdr_molecules): molecules is a
    list of (ID, Tmin, Tmax, n) in the order of molecules.input, n the number
    densities at the nspace depth points; nHmin the H- densities, the last
    record. Per molecule: counted ID (xdr short, xdr string); Nv, NJ; Ediss,
    Tmin, Tmax; n; if Nv > 0 and NJ > 0 the Nv x NJ energies. ValueError if
    the file does not parse."""
    b = read_bytes(path)
    mols = []
    try:
        nmol = struct.unpack(">i", b[:4])[0]
        if nmol <= 0:
            raise ValueError(f"Nmolecule {nmol}")
        off = 4
        for _ in range(nmol):
            ident, off = xdr_counted_string(b, off)
            if not ident:
                raise ValueError("empty molecule ID")
            nv, nj = struct.unpack(">2i", b[off:off + 8])
            tmin, tmax = struct.unpack(">2d", b[off + 16:off + 32])    # after Nv, NJ, Ediss
            off += 32
            n = struct.unpack(f">{nspace}d", b[off:off + 8 * nspace])
            off += 8 * nspace
            if nv > 0 and nj > 0:
                off += 8 * nv * nj
            mols.append((ident.strip(), tmin, tmax, n))
        nhmin = struct.unpack(f">{nspace}d", b[off:off + 8 * nspace])
        if off + 8 * nspace != len(b):
            raise ValueError(f"{len(b) - off - 8 * nspace} bytes left after nHmin")
    except (struct.error, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"{path}: cannot parse molecules.out ({e})") from None
    return mols, nhmin


def read_spectrum(path, nrays):
    """spectrum.out of an unpolarized run: Nspect (int), lambda, I
    [Nspect][nrays], big-endian XDR. Returns (lambda, rows of nrays values)."""
    raw = read_bytes(path)
    ns = struct.unpack(">i", raw[:4])[0]
    lam = struct.unpack(f">{ns}d", raw[4:4 + 8 * ns])
    pos = 4 + 8 * ns
    I = struct.unpack(f">{ns * nrays}d", raw[pos:pos + 8 * ns * nrays])
    return lam, [I[i * nrays:(i + 1) * nrays] for i in range(ns)]


def read_stokes_spectrum(path, nrays):
    """spectrum.out of a polarized run: Nspect (int), lambda, I [Nspect][nrays]
    (doubles), a 12-byte record boundary, Q, U, V. All big-endian XDR.
    Returns (lambda, [I, Q, U, V]) with each Stokes array as a list of rows
    [Nspect][nrays]."""
    raw = read_bytes(path)
    ns = struct.unpack(">i", raw[:4])[0]
    lam = struct.unpack(f">{ns}d", raw[4:4 + 8 * ns])
    size, pos = ns * nrays, 4 + 8 * ns
    I = struct.unpack(f">{size}d", raw[pos:pos + 8 * size])
    pos += 8 * size + 12
    QUV = [struct.unpack(f">{size}d", raw[pos + 8 * size * i:pos + 8 * size * (i + 1)])
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
    return struct.unpack(f">{m}d", raw[16:16 + 8 * m])


def read_solveray(path):
    """I, Q, U, V of a solveray output (spectrum_*): mux, mu (2 doubles), I,
    the int Nspect = 0 of ray.input, then Q, U, V, all big-endian XDR."""
    b = read_bytes(path)
    M = (len(b) - 20) // 32
    if 32 * M + 20 != len(b):
        raise ValueError(f"{path}: size {len(b)} does not fit the solveray output layout "
                         f"(32 M + 20 bytes)")
    off = 16 + 8 * M + 4
    return [struct.unpack(f">{M}d", b[o:o + 8 * M])
            for o in (16, off, off + 8 * M, off + 16 * M)]
