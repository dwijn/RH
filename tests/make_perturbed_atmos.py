#!/usr/bin/env python3
"""Make an inhomogeneous 2-D or 3-D RH XDR atmosphere from a 1-D column.

usage: make_perturbed_atmos.py IN_2D.atmos OUT.atmos NX [NY] [--dx KM]
                               [--amp A] [--vamp KMS]

The first column of the 2-D input (e.g. Atmos/FALC_03x82.atmos, which has a
chromosphere, transition region and 6 hydrogen populations) is replicated
over an NX (2-D) or NX x NY (3-D) grid and perturbed deterministically, so
that the result stresses NLTE convergence in multi-D geometry:

  * the chromosphere and transition region are shifted vertically by up to
    +-dz_max (scaled by AMP), so neighbouring columns have the transition
    region at different heights (strong horizontal gradients);
  * the chromospheric temperature is modulated by +-35% * AMP;
  * vertical velocities of up to +-VAMP km/s and horizontal ones of half that
    grow with height above the temperature minimum.

Horizontal structure is a sum of a few periodic modes with fixed phases, so
the output is reproducible and periodic in x (and y). Electron and hydrogen
densities follow the vertical shift but are not made consistent with the new
temperature; use SOLVE_NE / chemical equilibrium if that matters.
"""
import argparse
import struct
import sys

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("src"); p.add_argument("dst")
p.add_argument("nx", type=int); p.add_argument("ny", type=int, nargs="?")
p.add_argument("--dx", type=float, default=100.0, help="grid spacing [km]")
p.add_argument("--amp", type=float, default=1.0, help="perturbation scale")
p.add_argument("--vamp", type=float, default=8.0, help="max |vz| [km/s]")
a = p.parse_args()

b = open(a.src, "rb").read()
Nx0, Nz, NH = struct.unpack(">3i", b[:12])
bc = struct.unpack(">3i", b[12:24])
pos = 24 + 8 * Nx0
z = np.frombuffer(b, ">f8", Nz, pos).copy(); pos += 8 * Nz
col = []                    # 2-D arrays are stored with x fastest (k*Nx + l)
for _ in range(5):                                  # T, ne, vturb, vx, vz
    col.append(np.frombuffer(b, ">f8", Nx0 * Nz, pos).reshape(Nz, Nx0)[:, 0].copy())
    pos += 8 * Nx0 * Nz
nH = np.frombuffer(b, ">f8", NH * Nx0 * Nz, pos).reshape(NH, Nz, Nx0)[:, :, 0].copy()
pos += 8 * NH * Nx0 * Nz
if pos != len(b):
    sys.exit(f"unexpected {len(b) - pos} trailing bytes in {a.src}")
T0, ne0, vt0 = col[0], col[1], col[2]

three_d = a.ny is not None
nx, ny = a.nx, (a.ny if three_d else 1)
X, Y = np.meshgrid(np.arange(nx) / nx, np.arange(ny) / ny, indexing="ij")
tp = 2 * np.pi


def modes(seed):
    """Periodic horizontal pattern in [-1, 1], fixed by seed."""
    r = np.random.default_rng(seed)
    f = np.zeros_like(X)
    for kx, ky, w in [(1, 1, 1.0), (2, 1, 0.6), (1, 3, 0.4), (3, 2, 0.3)]:
        if not three_d:
            ky = 0
        f += w * np.cos(tp * (kx * X + ky * Y) + r.uniform(0, tp))
    return f / np.abs(f).max()


f_shift, f_temp, f_vz, f_vx, f_vy = (modes(s) for s in range(5))

ztmin = z[np.argmin(T0)]                           # temperature minimum height
ramp = np.clip((z - ztmin) / 1000.0, 0.0, 1.0)     # 0 below Tmin, 1 at +1000 km
dz_max = 250.0 * a.amp


def interp(q, zz, log=False):
    """Column quantity q evaluated at heights zz (z is decreasing)."""
    v = np.log(q) if log else q
    r = np.interp(zz, z[::-1], v[::-1])
    return np.exp(r) if log else r


shape = (nx, ny, Nz)
T = np.empty(shape); ne = np.empty(shape); vt = np.empty(shape)
nHp = np.empty((NH,) + shape)
for i in range(nx):
    for j in range(ny):
        zs = z - dz_max * f_shift[i, j] * ramp      # shifted sampling heights
        Tc = interp(T0, zs, log=True)
        chromo = (Tc < 2.0e4) * ramp
        T[i, j] = Tc * (1.0 + 0.35 * a.amp * f_temp[i, j] * chromo)
        ne[i, j] = interp(ne0, zs, log=True)
        vt[i, j] = interp(vt0, zs)
        for n in range(NH):
            nHp[n, i, j] = interp(np.maximum(nH[n], 1e-30), zs, log=True)
vz = a.vamp * f_vz[:, :, None] * ramp
vx = 0.5 * a.vamp * f_vx[:, :, None] * ramp
vy = 0.5 * a.vamp * f_vy[:, :, None] * ramp
vz = np.broadcast_to(vz, shape); vx = np.broadcast_to(vx, shape)
vy = np.broadcast_to(vy, shape)

be = lambda arr: np.ascontiguousarray(arr, ">f8").tobytes()
with open(a.dst, "wb") as f:
    if three_d:     # rhsc3d/readatmos_xdr.c: Nx Ny Nz NHydr, bc[2], dx dy z ...
        f.write(struct.pack(">4i", nx, ny, Nz, NH))
        f.write(struct.pack(">2i", *bc[1:]))
        f.write(struct.pack(">2d", a.dx, a.dx))
        f.write(be(z))
        # 3-D arrays are stored with x fastest: index = (k*Ny + j)*Nx + i
        tr = lambda q: q.transpose(2, 1, 0)
        for q in (T, ne, vt, vx, vy, vz):
            f.write(be(tr(q)))
        for n in range(NH):
            f.write(be(tr(nHp[n])))
    else:           # rhsc2d/readatmos_xdr.c: Nx Nz NHydr, bc[3], dx[Nx], z ...
        f.write(struct.pack(">3i", nx, Nz, NH))
        f.write(struct.pack(">3i", *bc))
        f.write(be(np.full(nx, a.dx)))
        f.write(be(z))
        # 2-D arrays are stored with x fastest: index = k*Nx + l
        for q in (T, ne, vt, vx, vz):
            f.write(be(q[:, 0].T))
        for n in range(NH):
            f.write(be(nHp[n, :, 0].T))
dims = f"{nx}x{ny}x{Nz}" if three_d else f"{nx}x{Nz}"
print(f"wrote {a.dst}: {dims}, T {T.min():.0f}-{T.max():.0f} K, "
      f"|vz| <= {np.abs(vz).max():.1f} km/s")
