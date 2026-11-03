# ISurface reference

Every surface implements `ISurface`. `S(u, v)` is the surface, `Su` is its derivative with respect to `u`, and so on.

## Properties

| Property | Meaning |
|---|---|
| `type` | Class name, such as `SphereSurface` |
| `domainU`, `domainV` | `{ start, end }` of each parameter |
| `degreeU`, `degreeV` | Polynomial degree, or `undefined` if there is none |
| `isClosedU`, `isClosedV` | The edges at the two ends of the domain coincide |
| `isPeriodicU`, `isPeriodicV` | The parameter wraps around |
| `enclosesVolume` | Each direction is closed, or has both end edges collapsed to points |
| `area` | Total surface area, cached |

A surface can be closed without being periodic, such as a `SweepSurface` of a circle whose ends happen to meet.

## Evaluation

| Method | Returns |
|---|---|
| `pointAt(u, v)` | Point |
| `derivativeAt(u, v, orderU = 1, orderV = 0)` | Partial derivative. Orders `(0, 0)` give the point. |
| `tangentUAt(u, v)`, `tangentVAt(u, v)` | Unit `Su` or `Sv` |
| `normalAt(u, v)` | Unit `Su x Sv` |
| `frameAt(u, v)` | `{ origin, tangentU, tangentV, normal }`. The tangents aren't generally perpendicular. |

At a joint, derivatives are right-handed in both parameters: the value on the side of the larger parameter.

### Singular points

Where `Su` or `Sv` vanishes, such as a sphere's pole or a cone's apex:

- `normalAt` returns the limit from inside the domain, computed from `Suv` where possible. A sphere's pole gives the true normal. A cone's apex has no single normal, so the result depends on `u`.
- `tangentUAt` and `tangentVAt` return the direction at a nearby regular point.
- Curvatures are those at a nearby regular point.

## Curvature

| Method | Returns |
|---|---|
| `firstFundamentalFormAt(u, v)` | `{ E, F, G }`: `Su.Su`, `Su.Sv`, `Sv.Sv` |
| `secondFundamentalFormAt(u, v)` | `{ L, M, N }`: `Suu.n`, `Suv.n`, `Svv.n` |
| `principalCurvaturesAt(u, v)` | `{ k1, k2, direction1, direction2 }`, with `k1 >= k2` |
| `gaussianCurvatureAt(u, v)` | `k1 k2` |
| `meanCurvatureAt(u, v)` | `(k1 + k2) / 2` |

A curvature is positive where the surface bends towards its normal, so reversing a surface negates both curvatures and the mean curvature, and keeps the Gaussian curvature. `direction2` is `normal x direction1`. Where the two curvatures are equal, every tangent is principal, and `direction1` is `tangentU`.

## Structure

| Method | Returns |
|---|---|
| `breakpointsU()`, `breakpointsV()` | Parameters where smoothness may be lost, including both ends |
| `isoCurveU(u)` | The curve at fixed `u`, running in `v`, as an `ICurve` |
| `isoCurveV(v)` | The curve at fixed `v`, running in `u` |
| `boundaryCurves()` | The four edges: `v` start, `u` end, `v` end, `u` start. Each runs in its own increasing parameter, so they don't form a loop. An edge collapsed to a point is a curve of zero length. |

## Measures

| Method | Returns |
|---|---|
| `areaBetween(u0, u1, v0, v1)` | Area of a parameter rectangle. The limits are clamped. |
| `enclosedVolume()` | Volume if `enclosesVolume`, otherwise `undefined`. It assumes the surface doesn't cross itself. |
| `centroid()` | Centroid of a thin sheet of uniform density |

## Sampling and queries

| Method | Returns |
|---|---|
| `sample(countU, countV)` | Grid of points, evenly spaced in the parameters, indexed `[u][v]` |
| `tessellate(options)` | `{ points, normals, parameters, triangles }`. Options: `chordTolerance` (default 0.01), `maxEdgeLength`, `maxPoints` (default 1 000 000). |
| `boundingBox()` | Axis-aligned box, cached |
| `closestPoint(point)` | `{ u, v, point, distance }` |
| `distanceTo(point)` | The distance alone |

## Intersections

| Method | Returns |
|---|---|
| `intersectCurve(curve, tolerance = 1e-6)` | `{ point, u, v, t, distance }[]`, sorted by `t`, the curve's parameter |
| `intersectLine(origin, direction, tolerance = 1e-6)` | The same, with `t` the signed distance from `origin` along `direction` |
| `intersectPlane(origin, normal, options)` | `{ points, parameters, closed }[]` |
| `intersect(other, options)` | `{ points, parameters, parametersOther, closed }[]` |

A closed polyline doesn't repeat its first point at the end. The options are `tolerance` (default `1e-6`) and `meshTolerance` (default 1/300 of the bounding-box diagonal; the smaller of the two for `intersect`). The [limits page](Algorithms.md#intersections) describes what these can miss.

## Shape

| Method | Returns |
|---|---|
| `isPlanar(tolerance = 1e-7)` | True if all sampled points lie in one plane |
| `planeNormal(tolerance = 1e-7)` | Unit normal of that plane on the side of `normalAt`, or `undefined` |

## Editing

| Method | Returns |
|---|---|
| `trim(u0, u1, v0, v1)` | The piece over a rectangle, with the **same parameters**. The range is clamped and must have a positive size. A periodic surface can't be trimmed across its seam. |
| `splitU(u)`, `splitV(v)` | The pieces before and after. The parameter must be strictly inside the domain. |
| `reversed()` | Same geometry with the normal flipped, by running `u` the other way: `reversed().pointAt(uStart + uEnd - u, v)` equals `pointAt(u, v)`. |
| `transformed(transform)` | The surface under an affine transform, with the **same parameters**: `transformed(T).pointAt(u, v)` equals `T(pointAt(u, v))` |
| `offset(distance)` | Parallel surface along the normal. Parameters are kept except on `ConeSurface`, where `v` shifts. |
| `toNurbs(tolerance = 1e-4)` | NURBS form. See the [catalogue](Surfaces.md#what-each-operation-returns) for which are exact. |

## Other

| Method | Returns |
|---|---|
| `equals(other, tolerance = 1e-7)` | True if both agree at a 9 by 9 grid of normalized parameters. It compares points, so the same geometry parameterized differently isn't equal. |
| `toJSON()` | A plain object for `surfaceFromJSON`. `ParametricSurface` throws. |

## Errors

Invalid input throws `RangeError`: a non-positive size, an empty domain, a `NaN` parameter, boundary curves that don't meet, and so on. `surfaceFromJSON` throws `TypeError` on an unknown type, and so does `toJSON` on a surface holding a function.
