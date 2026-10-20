#!/usr/bin/env python3
"""Cut a 2-D or 3-D RH XDR atmosphere out of a Bifrost snapshot.

usage: make_bifrost_atmos.py DIR OUT [--3d | --1d] [--y0 J] [--stride S]
                             [--x0 I] [--nx N] [--ny N]
                             [--zmin MM] [--zmax MM] [--zstride S]
                             [--vturb KMS]

DIR holds the public Bifrost FITS cubes of one snapshot (Carlsson et al.
2016, A&A 585, A4), e.g. en024048_hion: BIFROST_<run>_<var>_<snap>.fits
for lgtg, lgne, ux, uy, uz, bx, by, bz and the six hydrogen populations
lgn1..lgn6. These are in SI units, cell centred, with z the height (up) and
uz, bz positive upward; arrays are stored [z][y][x] with z increasing.

Writes OUT.atmos (the layout read by rhsc2d/ or rhsc3d/readatmos_xdr.c) and
OUT.B (B [T], gamma and chi [rad], read by readb_xdr.c). With --1d it writes
the single column (X0, Y0) as a text atmosphere on a height scale for rhf1d
(rhf1d/multiatmos.c), e.g. to compare a 2-D run with its columns in 1-D.

RH treats both horizontal boundaries as periodic (3-D: the only option;
2-D: FIXED would need the incoming intensity at every wavelength on the
sides). Bifrost is periodic in x and y, so

  * a full row (2-D) or the full domain (3-D) sampled every STRIDE points,
    with STRIDE a divisor of the box size (504 = 2^3 3^2 7), is periodic;
  * a window (--nx/--ny smaller than the strided box) is not: rays that
    leave one side re-enter at the other across a jump in the atmosphere.
    Fine for timing and memory, not for accuracy studies.

The 2-D case is the plane y = J, with x and z in the plane; By is kept, so
the field has its full 3-D orientation (chi measured from x towards y).
Heights are cut to [ZMIN, ZMAX] Mm and sampled every ZSTRIDE points; RH
wants z decreasing (top first). Bifrost has no microturbulence; VTURB
[km/s] is added uniformly.
"""
import argparse
import glob
import os
import re
import struct
import sys

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("dir"); p.add_argument("out")
p.add_argument("--3d", dest="three_d", action="store_true")
p.add_argument("--1d", dest="one_d", action="store_true",
               help="column X0 (strided) of row Y0, for rhf1d")
p.add_argument("--y0", type=int, default=0, help="row for the 2-D plane")
p.add_argument("--stride", type=int, default=1, help="horizontal sampling")
p.add_argument("--x0", type=int, default=0, help="window start (strided)")
p.add_argument("--y0w", type=int, default=0, help="3-D window start in y")
p.add_argument("--nx", type=int, help="window size (default: all)")
p.add_argument("--ny", type=int, help="3-D window size in y")
p.add_argument("--zmin", type=float, default=-0.5, help="bottom [Mm]")
p.add_argument("--zmax", type=float, default=3.0, help="top [Mm]")
p.add_argument("--zstride", type=int, default=1, help="vertical sampling")
p.add_argument("--vturb", type=float, default=0.0, help="[km/s]")
a = p.parse_args()


def fits_cube(path):
    """Memory map the primary array and read the 1-D z-coordinate extension
    (no astropy needed)."""
    with open(path, "rb") as f:
        cards = {}
        hdr = b""
        while b"END" not in cards:
            block = f.read(2880)
            hdr += block
            for i in range(0, 2880, 80):
                cards[block[i:i + 8].strip()] = block[i + 10:i + 80].split(b"/")[0].strip()
        if int(cards[b"BITPIX"]) != -32 or int(cards[b"NAXIS"]) != 3:
            sys.exit(f"{path}: expected a 3-D float32 cube")
        n1, n2, n3 = (int(cards[b"NAXIS%d" % i]) for i in (1, 2, 3))
        data = np.memmap(path, ">f4", "r", len(hdr), (n3, n2, n1))
        f.seek(len(hdr) + -(-4 * n1 * n2 * n3 // 2880) * 2880)
        ext = f.read()
    m = re.search(rb"END {77}", ext)
    start = -(-(m.end()) // 2880) * 2880
    z = np.frombuffer(ext, ">f4", n3, start).astype(float)
    return data, z, float(cards[b"CDELT1"]) * 1e3          # dx [km]


def path(var):
    hits = glob.glob(os.path.join(a.dir, f"BIFROST_*_{var}_*.fits"))
    if len(hits) != 1:
        sys.exit(f"need exactly one cube for {var} in {a.dir}, found {hits}")
    return hits[0]


T3, z, dx0 = fits_cube(path("lgtg"))
Nz0, Ny0, Nx0 = T3.shape
s = a.stride
if Nx0 % s or Ny0 % s:
    print(f"warning: stride {s} does not divide {Nx0}x{Ny0}; the strided "
          "grid is not periodic", file=sys.stderr)
kz = np.nonzero((z >= a.zmin) & (z <= a.zmax))[0][::-1][::a.zstride]  # top first
nxs = len(range(0, Nx0, s)); nys = len(range(0, Ny0, s))
nx = 1 if a.one_d else (a.nx or nxs)
ix = (np.arange(a.x0, a.x0 + nx) % nxs) * s
if a.three_d:
    ny = a.ny or nys
    iy = (np.arange(a.y0w, a.y0w + ny) % nys) * s
else:
    ny = 1
    iy = np.array([a.y0])
periodic = a.one_d or nx == nxs and ny == (nys if a.three_d else 1) and Nx0 % s == 0 \
    and (not a.three_d or Ny0 % s == 0)


def cut(var):
    """[z][y][x] subset, top first, as float64."""
    c = fits_cube(path(var))[0] if var != "lgtg" else T3
    return np.asarray(c[kz][:, iy][:, :, ix], float)


T = 10.0 ** cut("lgtg")
ne = 10.0 ** cut("lgne")
v = [cut(c) / 1e3 for c in ("ux", "uy", "uz")]                  # km/s
nH = [10.0 ** cut(f"lgn{n}") for n in range(1, 7)]
B = [cut(c) for c in ("bx", "by", "bz")]                        # T
Bmod = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
gamma = np.arccos(np.divide(B[2], Bmod, out=np.ones_like(Bmod), where=Bmod > 0))
chi = np.arctan2(B[1], B[0])
vt = np.full_like(T, a.vturb)
Nz = len(kz)
zkm = z[kz] * 1e3
dx = s * dx0

if a.one_d:     # rhf1d/multiatmos.c: text, height [km], ne and nH in cm^-3
    with open(a.out + ".atmos", "w") as f:
        f.write(f"* Bifrost column x = {ix[0]}, y = {iy[0]} (full-resolution indices)\n"
                f"  BIFROST_{ix[0]}_{iy[0]}\n  Height scale\n* log g [cm s^-2]\n  4.44\n"
                f"* Ndep\n  {Nz}\n*\n* height [km]   T [K]   ne [cm^-3]   V [km/s]   Vturb [km/s]\n")
        for k in range(Nz):
            f.write(f"  {zkm[k]:.8E}  {T[k, 0, 0]:.8E}  {ne[k, 0, 0] * 1e-6:.8E}"
                    f"  {v[2][k, 0, 0]:.8E}  {vt[k, 0, 0]:.8E}\n")
        f.write("*\n* Hydrogen populations [cm^-3]\n*  nh(1)  nh(2)  nh(3)  nh(4)  nh(5)  np\n")
        for k in range(Nz):
            f.write("  " + "  ".join(f"{q[k, 0, 0] * 1e-6:.8E}" for q in nH) + "\n")
    with open(a.out + ".B", "wb") as f:
        for q in (Bmod, gamma, chi):
            f.write(np.ascontiguousarray(q.ravel(), ">f8").tobytes())
    print(f"wrote {a.out}.atmos/.B: column of {Nz} points, z = {zkm[-1]:.0f}..{zkm[0]:.0f} km, "
          f"T = {T.min():.0f}..{T.max():.3g} K")
    sys.exit(0)

be = lambda q: np.ascontiguousarray(q, ">f8").tobytes()   # [z][y][x]: x fastest
bc_top, bc_bottom = 1, 2                                   # ZERO, THERMALIZED
with open(a.out + ".atmos", "wb") as f:
    if a.three_d:   # rhsc3d/readatmos_xdr.c: Nx Ny Nz NHydr, bc[2], dx dy z ...
        f.write(struct.pack(">4i", nx, ny, Nz, 6))
        f.write(struct.pack(">2i", bc_top, bc_bottom))
        f.write(struct.pack(">2d", dx, dx))
        f.write(be(zkm))
        for q in (T, ne, vt, v[0], v[1], v[2]):
            f.write(be(q))
    else:           # rhsc2d/readatmos_xdr.c: Nx Nz NHydr, bc[3], dx[Nx], z ...
        f.write(struct.pack(">3i", nx, Nz, 6))
        f.write(struct.pack(">3i", 1, bc_top, bc_bottom))   # PERIODIC
        f.write(be(np.full(nx, dx)))
        f.write(be(zkm))
        for q in (T, ne, vt, v[0], v[2]):
            f.write(be(q))
    for q in nH:
        f.write(be(q))
with open(a.out + ".B", "wb") as f:
    for q in (Bmod, gamma, chi):
        f.write(be(q))

dims = f"{nx}x{ny}x{Nz}" if a.three_d else f"{nx}x{Nz}"
print(f"wrote {a.out}.atmos/.B: {dims}, dx = {dx:.1f} km, "
      f"z = {zkm[-1]:.0f}..{zkm[0]:.0f} km, T = {T.min():.0f}..{T.max():.3g} K, "
      f"|B| <= {Bmod.max() * 1e4:.0f} G, "
      f"{'periodic' if periodic else 'NOT periodic (window)'}")
