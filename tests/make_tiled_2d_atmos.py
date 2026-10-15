#!/usr/bin/env python3
"""Tile the columns of a 2-D RH XDR atmosphere to make a larger grid.

usage: make_tiled_2d_atmos.py IN.atmos OUT.atmos NCOLS

Only used to get a grid large enough for memory-scaling checks from the
small atmosphere shipped in Atmos/. Handles the layout read by
rhsc2d/readatmos_xdr.c when no irradiated boundary data trails the file.
"""
import struct
import sys

import numpy as np

src, dst, ncols = sys.argv[1], sys.argv[2], int(sys.argv[3])
b = open(src, "rb").read()
Nx, Nz, NH = struct.unpack(">3i", b[:12])
bc = struct.unpack(">3i", b[12:24])
if bc[1] == 0 or bc[2] == 0:
    sys.exit("irradiated boundaries carry extra data that is not handled here")
pos = 24
dx = np.frombuffer(b, ">f8", Nx, pos); pos += 8 * Nx
z = np.frombuffer(b, ">f8", Nz, pos); pos += 8 * Nz
fields = []
for _ in range(5):                              # T, ne, vturb, vx, vz
    fields.append(np.frombuffer(b, ">f8", Nx * Nz, pos).reshape(Nx, Nz)); pos += 8 * Nx * Nz
nH = np.frombuffer(b, ">f8", NH * Nx * Nz, pos).reshape(NH, Nx, Nz); pos += 8 * NH * Nx * Nz
if pos != len(b):
    sys.exit(f"unexpected {len(b) - pos} trailing bytes in {src}")

reps = -(-ncols // Nx)
tile = lambda a: np.tile(a, (reps,) + (1,) * (a.ndim - 1))[:ncols]
with open(dst, "wb") as f:
    f.write(struct.pack(">3i", ncols, Nz, NH))
    f.write(struct.pack(">3i", *bc))
    f.write(np.asarray(np.resize(dx, ncols), ">f8").tobytes())
    f.write(np.asarray(z, ">f8").tobytes())
    for a in fields:
        f.write(np.asarray(tile(a), ">f8").tobytes())
    for n in range(NH):
        f.write(np.asarray(tile(nH[n]), ">f8").tobytes())
print(f"wrote {dst}: Nx={ncols} Nz={Nz} Nspace={ncols * Nz}")
