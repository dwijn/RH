"""Write a magnetic field for a 2-D or 3-D RH XDR atmosphere (STOKES_INPUT).

usage: make_field.py 2d|3d ATMOS OUT

The field varies across the columns so that Q and U change with the
reference frame from column to column: strength 500 G (1 + 0.5 sin phi),
inclination 45 deg + 30 deg sin phi and azimuth phi, with
phi = 2 pi l / Nx in 2-D and 2 pi (l / Nx + m / Ny) in 3-D, the same at
all heights. RH reads B in tesla, the angles in radians, as three XDR
double arrays of Nspace values (x fastest: [Nz][Nx] in 2-D,
[Nz][Ny][Nx] in 3-D), in the order B, inclination, azimuth (readb_xdr.c).
"""
import math
import struct
import sys

geom, src, dst = sys.argv[1:4]
with open(src, "rb") as f:
    head = f.read(16)
if geom == "2d":
    Nx, Nz = struct.unpack(">2i", head[:8])
    Ny = 1
elif geom == "3d":
    Nx, Ny, Nz = struct.unpack(">3i", head[:12])
else:
    sys.exit("geometry must be 2d or 3d")

phi = [2.0 * math.pi * (l / Nx + (m / Ny if geom == "3d" else 0.0))
       for m in range(Ny) for l in range(Nx)]
B = [0.05 * (1.0 + 0.5 * math.sin(p)) for p in phi] * Nz
gamma = [math.radians(45.0 + 30.0 * math.sin(p)) for p in phi] * Nz
chi = phi * Nz

with open(dst, "wb") as f:
    for a in (B, gamma, chi):
        f.write(struct.pack(">%dd" % len(a), *a))
print(f"wrote {dst}: Nx={Nx} Ny={Ny} Nz={Nz} Nspace={Nx * Ny * Nz}")
