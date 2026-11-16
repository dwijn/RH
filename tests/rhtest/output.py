"""Readers of RH's binary outputs (struct only). The formats live here."""
import struct


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def read_doubles(path):
    """A file of native-endian doubles (J.dat, Jgas.dat)."""
    raw = read_bytes(path)
    return struct.unpack(f"={len(raw) // 8}d", raw)


def read_pops(path):
    """NLTE populations of a pops.*.out file as a flat list (the layout
    analyze_run.read_pops detects): XDR, an unknown header, then Nlevel and
    Nspace as ints, the populations and the LTE populations as doubles."""
    b = read_bytes(path)
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return struct.unpack(f">{nl * ns}d", b[off + 8:off + 8 + 8 * nl * ns])
    raise ValueError(f"{path}: cannot parse as a pops file (no layout with Nlevel, Nspace "
                     f"and 2 x Nlevel x Nspace doubles fits the size {len(b)})")


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
