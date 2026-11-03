# Surface catalogue

Formulas below are in the surface's local `Frame3`. Ranges are given as `{ u: [start, end], v: [start, end] }`, and angles are in radians.

## Elementary surfaces

### `PlaneSurface`

`(u, v, 0)`, a rectangle of the frame's XY plane. `PlaneSurface.centred(origin, normal, sizeU, sizeV, xHint?)` builds one around a point. Infinite planes aren't represented.

### `CylinderSurface`

`(r cos u, r sin u, v)`. `u` defaults to a full turn, and `v` to `[0, 1]`. The normal points outwards.

### `ConeSurface`

`(R(v) cos u, R(v) sin u, v)`, with `R(v) = radius + v tan(semiAngle)`. The radius may reach zero, at the apex, but not go negative. `semiAngle` may be negative, for a cone narrowing as `v` grows.

### `SphereSurface`

`r (cos v cos u, cos v sin u, sin v)`: `u` is longitude and `v` latitude, so the poles are at `v = ±pi/2`. Defaults to the whole sphere. `SphereSurface.centredAt(centre, radius)` builds one aligned with the world axes.

### `TorusSurface`

`((R + r cos v) cos u, (R + r cos v) sin u, r sin v)`. Defaults to the whole torus, periodic both ways. Only the ring form is tested: `R > r` is expected, and a self-crossing torus has an unreliable `enclosedVolume`.

## Quadrics and screw surfaces

### `EllipsoidSurface`

`(a cos v cos u, b cos v sin u, c sin v)`.

### `EllipticParaboloidSurface`

`(a v cos u, b v sin u, c v^2)`, with `v >= 0` and the vertex at `v = 0`. The normal there points along `-Z`.

### `HyperbolicParaboloidSurface`

`(u, v, u^2/a^2 - v^2/b^2)`: the saddle. Its range must be given.

### `HyperboloidSurface`

One sheet: `(a cosh v cos u, b cosh v sin u, c sinh v)`. Two sheets (the upper one, with `sheets = 2`): `(a sinh v cos u, b sinh v sin u, c cosh v)`, with `v >= 0`.

### `HelicoidSurface`

`(v cos u, v sin u, pitch u / 2pi)`: a straight line turning about and climbing the axis, as on a screw thread. `u` may span several turns. Its range must be given.

## Free-form surfaces

### `NurbsSurface`

`new NurbsSurface(degreeU, degreeV, controlPoints, knotsU, knotsV, weights?)`, with control points and weights indexed `[u][v]`. Each direction needs `count + degree + 1` knots. The domain is `[knots[degree], knots[count]]`, so unclamped knots work. Reports `BSplineSurface` if all weights are equal.

Also offers `insertKnotU(u, times)`, `insertKnotV(v, times)` and `mapControlPoints(map)`. Its iso-curves are exact `NurbsCurve`s.

### `BSplineSurface`

A non-rational `NurbsSurface` with two extra constructors:

- `BSplineSurface.uniform(degreeU, degreeV, points)`: clamped, with even knots over `[0, 1]`
- `BSplineSurface.interpolate(points, degreeU = 3, degreeV = 3)`: passes through a grid of points, over `[0, 1]`

### `BezierSurface`

Any degrees, over `[0, 1]`. Pass `weights` for the rational form, reported as `RationalBezierSurface`.

## Patches

### `BilinearSurface`

The patch through four corners, over `[0, 1]`: `new BilinearSurface(p00, p10, p01, p11)`, named `p` + `u` + `v`. Twisted corners make it a saddle.

### `GridSurface`

Piecewise bilinear through a grid of points: the surface counterpart of a polyline. The parameters are the grid indices.

### `HermiteSurface`

Bicubic patches from a point, `Su`, `Sv` and optionally `Suv` at each grid point, with optional knots. C1 across patches. Without twists, these are Ferguson patches. Tangents are derivatives with respect to the parameter, so their size matters.

## Surfaces from curves

### `ExtrusionSurface`

`C(u) + v D`, with `v` in `[0, 1]` by default. Any affine transform keeps the class.

### `RevolutionSurface`

`new RevolutionSurface(curve, axisOrigin, axisDirection, uStart = 0, uEnd = 2pi)`. `u` is the right-handed angle about the axis, with `u = 0` at the profile as given, and `v` is the curve's parameter. A profile rising along the axis gets an outward normal.

### `RuledSurface`

`(1 - v) A(u) + v B(u)`, with both curves read over `[0, 1]`. Straight lines join points at the same normalized parameter, so the curves' parameterization shapes the surface.

### `LoftSurface`

Through a series of curves, each read over `u` in `[0, 1]`. Across them, `v` follows natural cubic splines of the curves, so it is C2 in `v`. `{ parameters }` sets where each curve lies in `v`. By default they are spaced by the average distance between curves, over `[0, 1]`. `{ closed: true }` continues smoothly back to the first curve, and needs at least three curves. Two curves give a ruled surface.

### `SweepSurface`

A profile carried along a path. `u` is the profile's parameter and `v` the path's. The profile is taken as positioned at the start of the path, so `pointAt(u, pathStart)` is the profile itself.

- `{ mode: 'rotation-minimizing' }`, the default, keeps the profile at a fixed attitude to the path, with as little twist as possible. Derivatives are numeric.
- `{ mode: 'translation' }` only moves the profile, so it never turns. Derivatives are exact.

On a closed non-planar path, a rotation-minimizing sweep doesn't generally return to its starting attitude, so the ends may not meet.

### `PipeSurface`

A tube of constant radius around a path: `P(v) + r (cos u N(v) + sin u B(v))`, with `(N, B)` rotation-minimizing and the normal outwards. It folds over itself where the radius exceeds the path's radius of curvature. Derivatives are numeric.

## Filling boundaries

### `CoonsSurface`

Bilinearly blended patch filling four curves, over `[0, 1]`: `new CoonsSurface(bottom, top, left, right)`. `bottom` and `top` run in `u`; `left` and `right` run in `v`. It matches the edges in position only, not in tangent. The constructor throws if the corners don't meet.

### `GordonSurface`

Through a network of crossing curves, over `[0, 1]`: `new GordonSurface(curvesU, curvesV, options)`. It is the sum of a loft through each family, minus the surface through their crossing points.

**Each `u` curve must cross each `v` curve at the same normalized parameters as its neighbours.** Curves that cross geometrically but with uneven parameters are refused with a `RangeError`, because the formula wouldn't then pass through them. `parametersU` and `parametersV` state the crossings; by default they are found from the first curve of each family.

## Compound and custom surfaces

### `CompositeSurface`

Patches joined into a grid, indexed `[u][v]`. Each keeps its parameter speed, so the domain is `[0, sum of patch sizes]` in each direction, and derivatives equal those of the patches. Patches in a column must have the same size in `u`, and those in a row the same size in `v`. The constructor throws on a gap.

### `ParametricSurface`

`new ParametricSurface(fn, uStart, uEnd, vStart, vEnd, options)`. Options: `derivative`, `breakpointsU`, `breakpointsV`, `periodicU` and `periodicV`. Without a derivative function, derivatives are numeric. It can't be saved to JSON.

## General forms

Operations return these when they can't keep the original class. Each wraps the original surface and keeps its parameters, except as noted.

| Class | Parameterization |
|---|---|
| `TrimmedSurface` | The same as the original, over a sub-rectangle |
| `ReversedSurface` | `u` mirrored over the same domain |
| `TransformedSurface` | The same as the original |
| `OffsetSurface` | The same as the original. It follows the base exactly, so it can fold over itself where the distance exceeds a radius of curvature. Derivatives are numeric. |

## What each operation returns

"Same" means the same class, in exact form.

| Surface | `trim` | `reversed` | `transformed` | `offset` | `toNurbs` |
|---|---|---|---|---|---|
| `PlaneSurface` | Same | `ReversedSurface` | Same, if similarity | Same | Exact, same parameters |
| `CylinderSurface`, `ConeSurface`, `SphereSurface`, `TorusSurface` | Same | `ReversedSurface` | Same, if similarity | Same, if the size stays positive | Exact |
| `EllipsoidSurface`, `EllipticParaboloidSurface`, `HyperboloidSurface` | Same | `ReversedSurface` | Same, if similarity | `OffsetSurface` | Exact |
| `HyperbolicParaboloidSurface` | Same | `ReversedSurface` | Same, if similarity | `OffsetSurface` | Exact, over `[0, 1]` |
| `HelicoidSurface` | Same | `ReversedSurface` | Same, if similarity | `OffsetSurface` | Approximate |
| `BezierSurface` | `NurbsSurface` | `NurbsSurface` | Same | `OffsetSurface` | Itself |
| `NurbsSurface`, `BSplineSurface` | `NurbsSurface` | `NurbsSurface` | `NurbsSurface` | `OffsetSurface` | Itself |
| `BilinearSurface`, `GridSurface` | `TrimmedSurface` | `ReversedSurface` | Same | `OffsetSurface` | Exact, same parameters |
| `HermiteSurface` | `TrimmedSurface` | `ReversedSurface` | Same | `OffsetSurface` | Exact, same parameters |
| `ExtrusionSurface` | `TrimmedSurface` | `ReversedSurface` | Same | `OffsetSurface` | As the curve's |
| `RevolutionSurface` | `TrimmedSurface` | `ReversedSurface` | Same, if rotation, translation or uniform scale | `OffsetSurface` | As the curve's |
| `RuledSurface`, `LoftSurface`, `CoonsSurface`, `GordonSurface`, `CompositeSurface` | `TrimmedSurface` | `ReversedSurface` | Same | `OffsetSurface` | Approximate |
| `SweepSurface`, `PipeSurface`, `ParametricSurface` | `TrimmedSurface` | `ReversedSurface` | `TransformedSurface` | `OffsetSurface` | Approximate |

A similarity transform is a rotation, reflection, translation or uniform scale. Other transforms return a `TransformedSurface` where the table says "if similarity", whose `toNurbs` is exact if the original's is.

"Exact" without "same parameters" means the same geometry at a different parameter: rational forms can't reproduce an angle parameter. "Approximate" means C1 piecewise bicubic with the same parameters. "As the curve's" means exact if the curve's `toNurbs` is exact, though the curve's parameter can change.

## What isn't included

- **Infinite planes, cylinders and cones**, because every surface must have a finite area. Use a large patch.
- **Trimmed surfaces with free-form boundaries**, as in a B-rep face. `trim` cuts a rectangle of parameters only. Curves on a surface and face topology belong in a solid-modelling layer.
- **Triangular Bezier patches and subdivision surfaces**, because `ISurface` maps from a rectangle. Use a `BezierSurface` with a collapsed edge, or a `GridSurface` of the subdivided mesh.
- **Rolling-ball fillets and blends between surfaces.** Those need a solved spine; a `CoonsSurface` over known edges covers simple cases.
- **Variable-radius pipes, and sweeps guided by a second rail.** Use a `LoftSurface` through placed sections, or a `ParametricSurface`.
- **NURBS degree elevation, knot removal, surface fitting and periodic NURBS as a distinct form.** Unclamped knots evaluate correctly, but `isPeriodicU` stays false, so parameters clamp.

For other parametric shapes, such as a catenoid, a Mobius strip or a superquadric, use `ParametricSurface`.
