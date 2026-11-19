# Stokes cone in the 2-D and 3-D formal solvers (B-63, N-47, Q-33)

Design note, 2026-11-18. Background: the interpolation audit (`/hao/solar1/dwijn/claude/RH/evidence/t16-interp-audit/REPORT.md`, items A1–A4) and the cone counts in `/hao/solar1/dwijn/claude/RH/evidence/t18-cone/`.

## Problem

Physical Stokes quantities lie in a cone:
- I ≥ |(Q, U, V)| for the Stokes vector and for the source vector S = ε/χ_I;
- χ_I ≥ |(χ_Q, χ_U, χ_V)|, i.e. |η′| ≤ 1, for the absorption part of K′.

The magneto-optical ρ′ has no bound.

The 2-D and 3-D short characteristics interpolate all of these at up- and downwind points, each component on its own. Where the weights are convex and shared by all components (bilinear, linear) the result stays in the cone: a convex combination of points in a convex cone. Where they are not, it can leave the cone:
- the 3-D XY bicubic has negative lobes and a clamp per component;
- the 2-D Bézier with Steffen-limited derivatives has weights that depend on each component.

Counted on modernize:
- **3-D B-63 reproducer:** |η′| up to 201 at 1.3×10⁻³ of the upwind points, the upwind I outside the cone at 2.7×10⁻³, I < 0 in the emergent H and K cores.
- **2-D slabs:** |η′| up to 2.4, the upwind I outside the cone up to 2.4 (3×10⁻⁴).
- **1-D:** none, since 1-D does not interpolate.

Independently, in all geometries the Bézier control points of S leave the cone, because each derivative is limited on its own (N-47; ≈10⁻⁵ of the control points in the runs, |Q|/I = 2.08 after one step in a constructed case).

## Design

**Interpolate in a space where the cone holds by construction.** Write each vector as a magnitude (I, S_I, or 1 for η′) times a polarization fraction p with |p| < 1, and map p to its rapidity

    w = atanh(|p|) p / |p|,        p = tanh(|w|) w / |w|.

This maps the open unit ball onto all of R³, so any interpolated w is a physical p. The name comes from the light-cone (Minkowski) form of the cone condition: χ = m (cosh φ, sinh φ n̂), with |η′| = tanh φ. With sin/cos one would either bound the angle again or fold back through |p| = 1, where 1 + K′ is singular.

- **Magnitude:** interpolated with the existing scalar scheme, which stays within the range of the nodes, so it is positive.
- **w and ρ′:** interpolated with weights that do not depend on the data: the unclamped bicubic in 3-D, and Bézier with the unlimited centred derivatives in 2-D. This is linear in the data, so the result turns with the (Q, U) frame. The frame dependence of N-26 goes away for these quantities.
- **Grid arrays:** made once per formal solution (wavelength and ray), with one atanh per grid point for S and for η′. Those of I are updated at each point as the sweep computes it, because 2-D stencils and the 3-D XZ/YZ planes read points of the current layer.
- **Where the old weights are already convex and shared** (3-D bilinear: the XZ/YZ planes always, every plane with `LINEAR_3D`; 2-D linear stencils), the component-wise interpolation is kept. It is in the cone, frame-independent and exact for linear fields, and those results do not change.

**Control points (N-47).** With both nodes in the cone, so is the secant control point a = S_node + (S_other − S_node)/3. A control point c outside the cone is moved along the segment from a to c until it reaches the cone (the first root of a quadratic), and the derivative is set from it.
- With K′ = 0 the DELO-Bézier step is a combination of the nodes and control points with weights ≥ 0, so it stays in the cone.
- The exponential step propagates a physical cubic S, so it stays in the cone too.
- Nothing changes where c is already in the cone, which includes any linear S.

## What it costs and what it gives up

- **Exactness for linear fields.** Rapidity interpolation is exact for a constant polarization fraction, not for a linear Stokes vector whose fraction varies. The 2-D linear-S check on meshes whose upwind points lie on columns is now 1.15×10⁻⁴, not 10⁻¹³. The convergence order for a cubic S is unchanged (2.8).
- **Cost:** about one atanh per grid point and formal solution, and one tanh per interpolated vector. K′ needs 3 (or 6 with magneto-optics) interpolations instead of 3 to 12.
- **Not covered:**
  - grid values already outside the cone (S_I < 0 of inverted or negative-opacity points, N-40) get w = 0 and keep their magnitude;
  - the DELO-Bézier step with K′ ≠ 0 is not a positive combination, so only the interpolated inputs and the control points are guaranteed.

## Stack

1. **Q-33** (results unchanged): `StokesGrid` passed down; one interpolation function per quantity and geometry.
2. **B-63:** rapidity interpolation in those functions, on the paths whose weights are not convex.
3. **N-47:** cone-aware control points in `Bezier3_Stokes_step`, all geometries.
