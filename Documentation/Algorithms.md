# Algorithms and limits

`SurfaceBase` implements `ISurface` numerically. Classes override it where a closed form exists. This page describes the numeric defaults and where they can fail.

Timings are from one desktop machine, on Node 26, and are only a guide.

## Derivatives

Analytic and NURBS surfaces have exact derivatives. So do surfaces built from curves, as far as the curves do. `SweepSurface` in rotation-minimizing mode, `PipeSurface`, `OffsetSurface` and `ParametricSurface` without a derivative function use finite differences: Fornberg weights on seven nodes per differentiated parameter, with a step of `1e-3` of the smooth span, shifted to stay inside it. Up to order 4 per parameter is available.

Rational surfaces differentiate their homogeneous form and apply the quotient rule (The NURBS Book, eq. 4.20).

## Normals and singular points

The normal is `Su x Sv`. A derivative counts as vanished below `1e-9` of its typical size over the surface, and the cross product as vanished below `1e-10` of `|Su| |Sv|`.

Near an edge where `Su` vanishes, `Su` is about `(v - v0) Suv`, so the limit of the normal comes from `Suv x Sv`, and likewise for `Sv`. If that also vanishes, the normal at a point `1e-4` of a span inside the domain is used.

**Limit:** this detects an edge collapsing to a point. A surface can also become singular at an isolated interior point, such as a badly folded offset, where the result is whatever the nudge finds.

Curvature at a singular point is measured at the same nudged point, because the forms themselves degenerate there.

## Area, volume and centroid

Each smooth cell between breakpoints is integrated by an 8 by 8 Gauss-Legendre rule, and compared with the same rule on its four quarters. A cell whose estimates differ by more than `1e-10` of the whole surface's estimate is divided, up to depth 6.

Volume is the flux `(1/3) S.(Su x Sv)`, by the divergence theorem, taken from the centre of the bounding box. **Limit:** it assumes that the surface doesn't cross itself, and that `enclosesVolume` means what it says: only the edges are checked, at nine samples.

`PlaneSurface`, `CylinderSurface`, `ConeSurface` and `SphereSurface` have closed-form areas.

## NURBS

Evaluation uses the basis-function derivatives of The NURBS Book, algorithm A2.3, in tensor product. `insertKnotU`, `insertKnotV` and `trim` apply the curve library's knot insertion (A5.1) to each line of control points, in homogeneous form, so rational shapes are kept exactly. `BSplineSurface.interpolate` follows A9.4, with chord-length parameters averaged over the lines of the grid. It solves dense systems, so cost grows up to cubically with each side of the grid; a 100 by 100 grid took about 50 ms.

The exact forms of the elementary surfaces and quadrics come from revolving a rational profile (A8.1), then applying an affine map where needed.

## Tessellation

The mesh is a tensor grid, so it has no cracks. Each smooth span starts as eight pieces. In each pass, a span of `u` or `v` is halved if any of its edges is further than `chordTolerance` from the surface at the 1/4, 1/2 or 3/4 point, or longer than `maxEdgeLength`. Cells also check the diagonal the triangles use, at its centre, which catches twist on a surface with straight edges. Passes stop at 16, or when the grid would exceed `maxPoints`.

Limits:

- **It checks samples**, so the true deviation can slightly exceed the tolerance. The tests allow 1.5 times the tolerance at points the algorithm didn't check. Use a smaller one if you need a margin.
- **Reaching `maxPoints` ends refinement quietly**, and the mesh may then miss the tolerance.
- **A local feature refines its whole row and column**, so one sharp detail on a large surface is expensive. A sphere of radius 2 gave about 4 000 points at a tolerance of `1e-2` in about 25 ms, and about 263 000 at `1e-4` in about 0.7 s.

## Bounding box

A grid of at least 24 intervals per direction, with at least 4 per smooth span. For each of the six directions, up to three local extremes of the grid are refined by alternating golden-section searches.

**Limit:** a feature narrower than the grid spacing can poke out of the box. The tests allow `2e-5` of the box's diagonal against a dense sample. Declare breakpoints on a `ParametricSurface` with such features.

## Closest point

A grid of at least 32 intervals per direction, with at least 4 per smooth span. Up to eight local minima are refined by damped Newton steps on the squared distance, falling back to Gauss-Newton where it isn't convex, with a parameter held at its limit when the step would push it out. The best result wins.

**Limit:** a minimum narrower than the grid spacing can be missed. `PlaneSurface` and `CylinderSurface` use closed forms, as does `SphereSurface` when the projected point lies within its domain, and don't have this limit.

## Intersections

All of these search a mesh and then refine, so what the mesh doesn't show, they can't find. Pass a smaller `meshTolerance` to look harder.

### Curves and lines

The curve is tessellated at the same tolerance as the surface, and a tree of bounding boxes finds nearby triangles. A segment piercing a triangle, or a segment's midpoint passing close over one, seeds a Levenberg-Marquardt refinement of `S(u, v) - C(t)`. A result counts if the final distance is within `tolerance`.

- **Tangent contacts aren't reliably found.** The close-pass seeds do catch clear cases, but the point is then accurate only to around the square root of the tolerance along the surface.
- **Close crossings merge.** Hits closer than `max(10 tolerance, 1e-4 of the smaller diagonal)` count as one.
- **Overlaps aren't segments.** A curve lying in the surface returns a scattering of points, and at most 2000 candidates are refined.

### Plane sections

The height above the plane is taken at each mesh point. Triangle edges that change sign are bisected on the true surface, and the resulting segments are chained into polylines. Heights within rounding noise of zero count as positive, so that sections along a seam aren't lost.

- **A closed curve smaller than the mesh can be missed**, as can a plane that only touches the surface.
- **A surface lying in the plane gives no section**, rather than its boundary.
- **Each point is exact, but the polyline's spacing follows the mesh**, not a chord tolerance of its own.

### Two surfaces

Triangle pairs from the two meshes give segments, which are chained into polylines. Each point is then pulled onto both surfaces by minimum-norm Gauss-Newton steps, and dropped if its gap stays above `tolerance`. Two spheres took about 45 ms, for a loop of about 100 points.

- **Tangent contacts aren't reported** where the meshes don't cross, and surfaces that touch along a curve can give fragments.
- **Overlaps aren't reported.** Coincident stretches give no result, or fragments at their edges, because coplanar triangles don't intersect.
- **Small loops and thin slivers smaller than the mesh can be missed.**
- **Branch points end a polyline**, so a figure of eight comes back as several pieces, and a section through a pole or seam may come back as several pieces rather than one loop.
- **`closed` describes the chain from the mesh.** If a point was dropped during refinement, the loop has a small gap where it was.
- **Parameters on a seam or pole can come from either side.**

## NURBS approximation

For classes without an exact form, `toNurbs` builds C1 piecewise bicubics from the point, tangents and twist on a grid, doubling the pieces per smooth span until the deviation is within tolerance.

**Limit:** deviation is measured at nine points per cell, and refinement stops at about 128 cells per direction, when the result is returned even if it misses the tolerance. Compare the result against the surface if this matters, for instance at a very small tolerance on a large surface. Cost is high for surfaces with numeric derivatives: a pipe around a helix took about 0.3 s at `1e-3`, and about 4 s at `1e-5`.

## Offsets

`PlaneSurface`, `CylinderSurface`, `ConeSurface`, `SphereSurface` and `TorusSurface` have exact offsets while their size stays positive. `OffsetSurface` follows `S + d n` directly. It doesn't remove folds where the distance exceeds a radius of curvature, and doesn't fill gaps at creases such as a `CompositeSurface`'s joints.

## Sweeps and pipes

Frames come from the double-reflection method of Wang et al. (2008), from a table of about 256 entries along the path, with each query stepping from its nearest entry.

**Limit:** around a closed non-planar path, frames don't return to their start, so the ends of a sweep or pipe may not meet, and may twist. No correction is applied. At a corner of the path, the profile turns abruptly, leaving a gap or overlap on the surface.

## Numerical tolerances

| Value | Use |
|---|---|
| `1e-9`, plus `1e-12` of the point's distance from the origin | Edges coincide (`isClosedU`, `isClosedV`, `enclosesVolume`) |
| `1e-7` | Planarity; `equals`; corners of a `CoonsSurface`; joints of a `CompositeSurface`, relative to its patches' size |
| `1e-6` | Default intersection tolerance; crossings of a `GordonSurface`, relative to its size |
| `1e-8` of the diagonal | Merging end points when chaining intersection polylines |

The tests run every surface at model sizes of `1e-3` and `1e6`, and areas, boxes, curvatures, closest points and tessellation hold up at both. Pass tolerances in proportion to the model where a method accepts them, particularly the intersection tolerance and `chordTolerance`, and rescale models well outside that range.
