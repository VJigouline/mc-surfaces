# MeerCADSurfaces

CAD surfaces in TypeScript. Planes, quadrics, NURBS, and swept, lofted and blended surfaces all implement one `ISurface` interface, so the same code can evaluate, measure, tessellate, intersect, trim, transform and save any of them.

- Works with [MeerCADCurves](https://github.com/VJigouline/MeerCADCurves): any of its curves can be extruded, revolved, lofted or swept
- Immutable surfaces: operations return new surfaces
- Closed-form results where they exist, and tested numerics where they don't
- Ships as an ES module with type declarations

## Install

The package isn't published to npm yet. Install it from GitHub, which also fetches and builds its `meercad-curves` dependency:

```sh
npm install github:VJigouline/MeerCADSurfaces
```

Node.js 22 or later is required.

## Quick example

```ts
import { CircularArc, Line } from 'meercad-curves';
import {
  ExtrusionSurface,
  Frame3,
  LoftSurface,
  SphereSurface,
  Transform3,
  Vector3,
  surfaceFromJSON,
} from 'meercad-surfaces';

// Every surface shares the ISurface interface.
const sphere = SphereSurface.centredAt(new Vector3(0, 0, 0), 2);
const tube = new ExtrusionSurface(CircularArc.circle(new Vector3(0, 0, -5), 1), new Vector3(0, 0, 10));

sphere.area;                                // 50.265...
sphere.enclosedVolume();                    // 33.510...
sphere.normalAt(0.5, 0.2);                  // unit Vector3
sphere.gaussianCurvatureAt(0.5, 0.2);       // 0.25
tube.tessellate({ chordTolerance: 0.01 });  // points, normals, parameters, triangles

// Queries
sphere.closestPoint(new Vector3(5, 0, 0)).distance;       // 3
sphere.intersectLine(new Vector3(-5, 0, 0), Vector3.X);   // two hits
sphere.intersect(tube);                                   // two closed loops

// Surfaces from curves
const loft = new LoftSurface([
  new Line(new Vector3(0, 0, 0), new Vector3(4, 0, 0)),
  new Line(new Vector3(0, 2, 1), new Vector3(4, 2, 2)),
  new Line(new Vector3(0, 4, 0), new Vector3(4, 4, 0)),
]);

// Operations return new surfaces, and keep the parameters.
const [left, right] = loft.splitU(0.4);
const moved = sphere.transformed(Transform3.translate(new Vector3(10, 0, 0)));
const nurbs = sphere.toNurbs();             // exact rational form

// Save and restore
const copy = surfaceFromJSON(JSON.parse(JSON.stringify(moved.toJSON())));
```

## Surfaces

| Group | Classes |
|---|---|
| Elementary | `PlaneSurface`, `CylinderSurface`, `ConeSurface`, `SphereSurface`, `TorusSurface` |
| Quadrics and screw | `EllipsoidSurface`, `EllipticParaboloidSurface`, `HyperbolicParaboloidSurface`, `HyperboloidSurface`, `HelicoidSurface` |
| Free-form | `BezierSurface` (with rational form), `BSplineSurface`, `NurbsSurface` |
| Patches | `BilinearSurface`, `GridSurface`, `HermiteSurface` |
| From curves | `ExtrusionSurface`, `RevolutionSurface`, `RuledSurface`, `LoftSurface`, `SweepSurface`, `PipeSurface` |
| Filling boundaries | `CoonsSurface`, `GordonSurface` |
| Compound | `CompositeSurface` |
| Custom | `ParametricSurface`, for any function of two parameters |

Operations that can't keep a surface's class return one of these general forms: `TrimmedSurface`, `ReversedSurface`, `TransformedSurface` or `OffsetSurface`.

## Documentation

- [Getting started](Documentation/GettingStarted.md): concepts and common tasks
- [ISurface reference](Documentation/ISurface.md): every property and method
- [Surface catalogue](Documentation/Surfaces.md): each surface's formula, parameters and exact operations
- [Algorithms and limits](Documentation/Algorithms.md): how results are computed, and where they stop being reliable

## Development

| Command | Purpose |
|---|---|
| `npm test` | Run the tests |
| `npm run coverage` | Run the tests with a coverage report |
| `npm run typecheck` | Type-check the sources and tests |
| `npm run build` | Compile to `dist` |

The core of the test suite is `tests/contract.test.ts`, which checks every surface in `tests/catalog.ts` against the whole interface, and `tests/scale.test.ts` repeats key checks on very small and very large models. Add a new surface to the catalogue and it is held to the same contract. `tests/surfaces.test.ts` compares results with known closed forms.

## License

[MIT](LICENSE)
