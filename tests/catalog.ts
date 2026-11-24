import {
  BezierCurve,
  BSplineCurve,
  CircularArc,
  CubicSpline,
  EllipticalArc,
  Helix,
  Line,
  Polyline,
  type ICurve,
} from 'mc-curves';
import {
  BezierSurface,
  BilinearSurface,
  BSplineSurface,
  CompositeSurface,
  ConeSurface,
  CoonsSurface,
  CylinderSurface,
  EllipsoidSurface,
  EllipticParaboloidSurface,
  ExtrusionSurface,
  Frame3,
  GordonSurface,
  GridSurface,
  HelicoidSurface,
  HermiteSurface,
  HyperbolicParaboloidSurface,
  HyperboloidSurface,
  LoftSurface,
  NurbsSurface,
  ParametricSurface,
  PipeSurface,
  PlaneSurface,
  RevolutionSurface,
  RuledSurface,
  SphereSurface,
  SweepSurface,
  TorusSurface,
  Transform3,
  Vector3,
  type ISurface,
} from '../src/index.js';
import { v } from './helpers.js';

/** A frame that is not aligned with the world, so tests catch mixed-up local and world coordinates. */
export const tilted = Frame3.fromNormal(v(1, 2, 3), v(1, 1, 2), v(1, -1, 0));

const wave = (x: number, y: number) => Math.sin(x) * Math.cos(y);
const net = (nx: number, ny: number) =>
  Array.from({ length: nx }, (_, i) => Array.from({ length: ny }, (_, j) => v(i, j, wave(i * 1.3, j * 0.9))));

const profile = new BezierCurve([v(2, 0, 0), v(3, 0, 1), v(1.5, 0, 2), v(2, 0, 3)]);
const wall = Frame3.fromAxes(v(0, 1, 0), Vector3.X, Vector3.Z, Vector3.Y.negate());
const plateA = new PlaneSurface(Frame3.WORLD, { u: [0, 1], v: [0, 1] });
const plateB = new PlaneSurface(Frame3.WORLD, { u: [1, 2], v: [0, 1] });
const plateC = new PlaneSurface(wall, { u: [0, 1], v: [0, 1] });
const plateD = new PlaneSurface(wall, { u: [1, 2], v: [0, 1] });

const networkU: ICurve[] = [0, 1, 2].map((j) => new BezierCurve([v(0, j, 0), v(1, j, j), v(2, j, 0)]));
const networkV: ICurve[] = [0, 0.5, 1].map((s) => {
  const lift = (j: number) => 2 * s * (1 - s) * j;
  return new BezierCurve([v(2 * s, 0, lift(0)), v(2 * s, 1, lift(1)), v(2 * s, 2, lift(2))]);
});

/** One instance of every surface kind, used by the contract tests. */
export const catalog: Record<string, () => ISurface> = {
  PlaneSurface: () => new PlaneSurface(tilted, { u: [-2, 3], v: [1, 4] }),
  CylinderSurface: () => new CylinderSurface(tilted, 2, { v: [-1, 3] }),
  'CylinderSurface (patch)': () => new CylinderSurface(tilted, 2, { u: [0.5, 2.5], v: [0, 2] }),
  ConeSurface: () => new ConeSurface(tilted, 1, 0.4, { v: [0, 3] }),
  'ConeSurface (apex)': () => new ConeSurface(tilted, 0, 0.5, { v: [0, 2] }),
  SphereSurface: () => new SphereSurface(tilted, 2.5),
  'SphereSurface (patch)': () => new SphereSurface(tilted, 2.5, { u: [0.3, 2], v: [-0.5, 1] }),
  TorusSurface: () => new TorusSurface(tilted, 3, 1),
  'TorusSurface (patch)': () => new TorusSurface(tilted, 3, 1, { u: [0, 2], v: [1, 4] }),
  EllipsoidSurface: () => new EllipsoidSurface(tilted, 3, 2, 1),
  EllipticParaboloidSurface: () => new EllipticParaboloidSurface(tilted, 2, 1, 0.5, { v: [0, 2] }),
  HyperbolicParaboloidSurface: () => new HyperbolicParaboloidSurface(tilted, 2, 3, { u: [-2, 2], v: [-1.5, 2] }),
  'HyperboloidSurface (one sheet)': () => new HyperboloidSurface(tilted, 2, 1, 1.5),
  'HyperboloidSurface (two sheets)': () => new HyperboloidSurface(tilted, 2, 1, 1.5, { v: [0, 1.2] }, 2),
  HelicoidSurface: () => new HelicoidSurface(tilted, 1.5, { u: [0, 9], v: [0.5, 2] }),
  BezierSurface: () => new BezierSurface(net(4, 3)),
  RationalBezierSurface: () =>
    new BezierSurface(net(3, 3), [
      [1, 2, 1],
      [0.5, 3, 1],
      [1, 1, 2],
    ]),
  BSplineSurface: () => BSplineSurface.uniform(3, 2, net(6, 5)),
  'BSplineSurface (interpolated)': () => BSplineSurface.interpolate(net(5, 5)),
  NurbsSurface: () =>
    new NurbsSurface(
      2,
      2,
      net(4, 4),
      [0, 0, 0, 1, 3, 3, 3],
      [0, 0, 0, 2, 3, 3, 3],
      Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) => 1 + ((i * j) % 3))),
    ),
  BilinearSurface: () => new BilinearSurface(v(0, 0), v(3, 0, 1), v(0, 2, 1), v(3, 2)),
  GridSurface: () => new GridSurface(net(4, 4)),
  HermiteSurface: () =>
    new HermiteSurface({
      points: net(3, 3),
      tangentsU: net(3, 3).map((row) => row.map((p) => v(1, 0, Math.cos(p.x * 1.3)))),
      tangentsV: net(3, 3).map((row) => row.map((p) => v(0, 1, Math.sin(p.y * 0.9)))),
      twists: net(3, 3).map((row) => row.map((p) => v(0, 0, 0.2 * p.x - 0.1 * p.y))),
      knotsU: [0, 1, 2.5],
      knotsV: [0, 2, 3],
    }),
  ExtrusionSurface: () => new ExtrusionSurface(BSplineCurve.uniform(3, [v(0, 0), v(1, 3), v(3, -2, 1), v(5, 3), v(7, 0)]), v(0.5, 1, 4)),
  'ExtrusionSurface (of a circle)': () => new ExtrusionSurface(CircularArc.circle(v(1, 1), 2), v(0, 0, 3)),
  'ExtrusionSurface (of a polyline)': () => new ExtrusionSurface(new Polyline([v(0, 0), v(2, 0), v(2, 2), v(4, 3)]), v(0, 0, 2)),
  RevolutionSurface: () => new RevolutionSurface(profile, v(0, 0, 0), v(0, 0, 1)),
  'RevolutionSurface (part)': () => new RevolutionSurface(new Line(v(3, 0, 0), v(1, 1, 2)), v(1, 0, 0), v(0.3, 0, 1), 0.5, 3),
  RuledSurface: () => new RuledSurface(new EllipticalArc(tilted, 3, 2, 0, 2), new Line(v(0, 0, 5), v(4, 1, 6))),
  LoftSurface: () =>
    new LoftSurface([
      new BezierCurve([v(0, 0, 0), v(1, 1, 0), v(2, 0, 0)]),
      new BezierCurve([v(0, 0, 1), v(1, 2, 1.5), v(2, 0, 1)]),
      new BezierCurve([v(0, 0, 2.5), v(1, 0.5, 2), v(2, 0, 2.5)]),
      new BezierCurve([v(0, 0, 4), v(1, 1, 4), v(2, 0, 4)]),
    ]),
  'LoftSurface (closed)': () =>
    new LoftSurface(
      [0, 1, 2, 3].map((k) => {
        const a = (k * Math.PI) / 2;
        const r = k % 2 === 0 ? 2 : 3;
        return new Line(v(r * Math.cos(a), r * Math.sin(a), 0), v(r * Math.cos(a), r * Math.sin(a), 2));
      }),
      { closed: true },
    ),
  CoonsSurface: () =>
    new CoonsSurface(
      new BezierCurve([v(0, 0, 0), v(1, -0.5, 1), v(2, 0, 0)]),
      new BezierCurve([v(0, 2, 0), v(1, 2.5, -1), v(2, 2, 0)]),
      new BezierCurve([v(0, 0, 0), v(-0.5, 1, 0.5), v(0, 2, 0)]),
      new BezierCurve([v(2, 0, 0), v(2.5, 1, -0.5), v(2, 2, 0)]),
    ),
  GordonSurface: () => new GordonSurface(networkU, networkV),
  SweepSurface: () => new SweepSurface(new BezierCurve([v(0, 0, 0), v(1, 0, 1), v(2, 0, 0)]), new BezierCurve([v(0, 0, 0), v(0, 3, 0), v(3, 5, 1), v(5, 5, 3)])),
  'SweepSurface (translation)': () =>
    new SweepSurface(new Line(v(0, 0, 0), v(1, 0, 0)), new BezierCurve([v(0, 0, 0), v(0, 2, 1), v(1, 3, 3)]), { mode: 'translation' }),
  PipeSurface: () => new PipeSurface(new Helix(tilted, 4, 3, 1.5), 0.5),
  'PipeSurface (spline)': () => new PipeSurface(new CubicSpline([v(0, 0), v(5, 2), v(8, 0, 3), v(12, 4, 2)]), 0.4),
  CompositeSurface: () =>
    new CompositeSurface([
      [plateA, plateC],
      [plateB, plateD],
    ]),
  ParametricSurface: () => new ParametricSurface((a, b) => v(a, b, Math.sin(a) * Math.cos(b)), 0, 3, -1, 2),
  'TrimmedSurface (of a loft)': () => catalog['LoftSurface']!().trim(0.2, 0.8, 0.1, 0.7),
  'ReversedSurface (of a torus)': () => catalog['TorusSurface']!().reversed(),
  'TransformedSurface (stretched sphere)': () => catalog['SphereSurface']!().transformed(Transform3.scale(1, 2, 0.5)),
  'OffsetSurface (of a paraboloid)': () => new EllipticParaboloidSurface(tilted, 2, 1, 0.5, { v: [0.3, 1.5] }).offset(0.2),
};

/** Surfaces that hold a function, so they can't be saved to JSON. */
export const notSerializable = new Set(['ParametricSurface']);
