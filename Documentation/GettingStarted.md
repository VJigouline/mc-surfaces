# Getting started

## Concepts

A surface is a map from a rectangle of parameters `(u, v)` to 3D space. `domainU` and `domainV` give the rectangle. Infinite planes and cylinders aren't represented, because every surface needs a finite area.

The math types (`Vector3`, `Frame3`, `Transform3`, `BoundingBox3`) come from `meercad-curves`, and this package re-exports them, so both libraries share one set. Curves come from `meercad-curves` directly.

Three rules apply to every surface:

- **Surfaces are immutable.** `trim`, `transformed` and the like return new surfaces.
- **Out-of-range parameters don't throw.** They are clamped to the domain, except in a periodic direction, where they wrap around.
- **The normal is `Su x Sv`.** So it depends on how the surface is parameterized, and a mirror turns it inside out. Use `reversed()` to flip it.

Analytic surfaces live in a local `Frame3`. A surface built from a tilted frame is tilted:

```ts
import { CylinderSurface, Frame3, Vector3 } from 'meercad-surfaces';

const frame = Frame3.fromNormal(new Vector3(1, 2, 3), new Vector3(0, 1, 1));
const cylinder = new CylinderSurface(frame, 2, { v: [0, 10] });
```

## Evaluating

```ts
cylinder.pointAt(1, 5);
cylinder.derivativeAt(1, 5);        // Su
cylinder.derivativeAt(1, 5, 0, 1);  // Sv
cylinder.derivativeAt(1, 5, 1, 1);  // Suv
cylinder.normalAt(1, 5);
cylinder.frameAt(1, 5);             // origin, tangentU, tangentV, normal
```

## Curvature

```ts
const { k1, k2, direction1, direction2 } = cylinder.principalCurvaturesAt(1, 5);
cylinder.gaussianCurvatureAt(1, 5);
cylinder.meanCurvatureAt(1, 5);
```

A curvature is positive where the surface bends towards its normal. A sphere's normal points outwards, so its curvatures are `-1/r`.

## Building surfaces from curves

```ts
import { BezierCurve, CircularArc, Line } from 'meercad-curves';
import { LoftSurface, PipeSurface, RevolutionSurface, RuledSurface, SweepSurface, Vector3 } from 'meercad-surfaces';

const profile = new BezierCurve([new Vector3(2, 0, 0), new Vector3(3, 0, 1), new Vector3(2, 0, 2)]);
const vase = new RevolutionSurface(profile, Vector3.ZERO, Vector3.Z);

const path = new BezierCurve([new Vector3(0, 0, 0), new Vector3(0, 5, 0), new Vector3(5, 5, 2)]);
const hose = new PipeSurface(path, 0.25);
```

Curves given to `RuledSurface`, `LoftSurface`, `CoonsSurface` and `GordonSurface` are read over `[0, 1]`, whatever their own domains, and must run the same way. Mismatched directions twist the surface.

## Tessellating for display

```ts
const mesh = vase.tessellate({ chordTolerance: 0.01 });
mesh.points;      // Vector3[]
mesh.normals;     // Vector3[], one per point
mesh.parameters;  // { u, v }[], one per point
mesh.triangles;   // three indices per triangle
```

Points aren't shared across a seam or at a pole, so that each keeps its own parameters. Triangle winding agrees with the normal.

## Queries

```ts
vase.closestPoint(new Vector3(5, 0, 1));           // { u, v, point, distance }
vase.boundingBox();
vase.intersectLine(new Vector3(9, 0, 1), Vector3.X); // { point, u, v, t, distance }[]
vase.intersectCurve(path);
vase.intersectPlane(new Vector3(0, 0, 1), Vector3.Z); // polylines
vase.intersect(hose);                                 // polylines with parameters on both surfaces
```

Read [the limits of each](Algorithms.md) before relying on them. In particular, tangent contacts are not reliably found.

## Measuring

```ts
vase.area;
vase.areaBetween(0, 1, 0, 0.5);
vase.centroid();
vase.enclosesVolume;     // needs both directions closed or collapsed to a point
vase.enclosedVolume();   // undefined unless enclosesVolume is true
```

A vase with open ends doesn't enclose a volume, so `enclosedVolume()` is `undefined`. A sphere, torus or full ellipsoid does.

## Editing

```ts
const part = vase.trim(0, 3, 0.2, 0.8);    // keeps the parameters
const [a, b] = vase.splitV(0.5);
const flipped = vase.reversed();           // same geometry, opposite normal
const moved = vase.transformed(transform); // keeps the parameters
const shell = vase.offset(0.1);
const nurbs = vase.toNurbs(1e-4);
```

## Writing your own surface

The quickest way is `ParametricSurface`:

```ts
import { ParametricSurface, Vector3 } from 'meercad-surfaces';

const wave = new ParametricSurface((u, v) => new Vector3(u, v, Math.sin(u) * Math.cos(v)), 0, 6, 0, 6);
```

Derivatives are numeric unless you pass a `derivative` function. Declare `breakpointsU` and `breakpointsV` wherever it isn't smooth, and don't expect it to save to JSON.

For a class of your own, extend `SurfaceBase` and implement `domainU`, `domainV`, `evaluate` and `toJSON`. Everything else in `ISurface` then works numerically, and you can override whatever has a closed form. To reuse the tests, add it to `tests/catalog.ts`.

## Saving

```ts
import { surfaceFromJSON } from 'meercad-surfaces';

const text = JSON.stringify(vase.toJSON());
const copy = surfaceFromJSON(JSON.parse(text));
```

A surface built on a `ParametricCurve`, or a `ParametricSurface`, can't be saved, because each holds a function. `toJSON` throws for those. To read an iso-curve back, use this package's `curveFromJSON`, which adds that one type to the curve library's reader; it doesn't reach an iso-curve inside a curve wrapper such as a `TrimmedCurve`.
