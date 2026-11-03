import { curveFromJSON as readPlainCurve, Frame3, Transform3, Vector3, type CurveJSON, type ICurve } from 'meercad-curves';
import { CompositeSurface } from './Compound.js';
import { CoonsSurface, ExtrusionSurface, GordonSurface, LoftSurface, RevolutionSurface, RuledSurface } from './CurveSurfaces.js';
import { ConeSurface, CylinderSurface, PlaneSurface, SphereSurface, TorusSurface } from './Elementary.js';
import type { ISurface, SurfaceJSON } from './ISurface.js';
import { BezierSurface, BSplineSurface, NurbsSurface } from './NurbsSurface.js';
import { BilinearSurface, GridSurface, HermiteSurface } from './Patches.js';
import { EllipsoidSurface, EllipticParaboloidSurface, HelicoidSurface, HyperbolicParaboloidSurface, HyperboloidSurface } from './Quadrics.js';
import { IsoCurve, OffsetSurface, ReversedSurface, TransformedSurface, TrimmedSurface } from './SurfaceBase.js';
import { PipeSurface, SweepSurface } from './SweepSurfaces.js';

type Json = Record<string, any>;

const point = (a: readonly number[]) => Vector3.fromArray(a);
const grid = (a: readonly (readonly (readonly number[])[])[]) => a.map((row) => row.map(point));

function frameOf(j: Json): Frame3 {
  const f = j['frame'];
  return Frame3.fromAxes(point(f.origin), point(f.x), point(f.y), point(f.z));
}

/** Like the curve library's reader, but it also knows the iso-curves that surfaces produce. */
export function curveFromJSON(json: CurveJSON): ICurve {
  if (json.type === 'IsoCurve') {
    const j = json as Json;
    return new IsoCurve(surfaceFromJSON(j['surface']), j['fixed'], j['value']);
  }
  return readPlainCurve(json);
}

const readers: Record<string, (j: Json) => ISurface> = {
  PlaneSurface: (j) => new PlaneSurface(frameOf(j), { u: j['u'], v: j['v'] }, j['scale']),
  CylinderSurface: (j) => new CylinderSurface(frameOf(j), j['radius'], { u: j['u'], v: j['v'] }, j['scale']),
  ConeSurface: (j) => new ConeSurface(frameOf(j), j['radius'], j['semiAngle'], { u: j['u'], v: j['v'] }, j['scale']),
  SphereSurface: (j) => new SphereSurface(frameOf(j), j['radius'], { u: j['u'], v: j['v'] }),
  TorusSurface: (j) => new TorusSurface(frameOf(j), j['majorRadius'], j['minorRadius'], { u: j['u'], v: j['v'] }),
  EllipsoidSurface: (j) => new EllipsoidSurface(frameOf(j), j['radiusX'], j['radiusY'], j['radiusZ'], { u: j['u'], v: j['v'] }),
  EllipticParaboloidSurface: (j) => new EllipticParaboloidSurface(frameOf(j), j['a'], j['b'], j['c'], { u: j['u'], v: j['v'] }),
  HyperbolicParaboloidSurface: (j) => new HyperbolicParaboloidSurface(frameOf(j), j['a'], j['b'], { u: j['u'], v: j['v'] }, j['scale']),
  HyperboloidSurface: (j) => new HyperboloidSurface(frameOf(j), j['a'], j['b'], j['c'], { u: j['u'], v: j['v'] }, j['sheets']),
  HelicoidSurface: (j) => new HelicoidSurface(frameOf(j), j['pitch'], { u: j['u'], v: j['v'] }, j['scale']),
  BezierSurface: (j) => new BezierSurface(grid(j['controlPoints'])),
  RationalBezierSurface: (j) => new BezierSurface(grid(j['controlPoints']), j['weights']),
  BSplineSurface: (j) => new BSplineSurface(j['degreeU'], j['degreeV'], grid(j['controlPoints']), j['knotsU'], j['knotsV']),
  NurbsSurface: (j) => new NurbsSurface(j['degreeU'], j['degreeV'], grid(j['controlPoints']), j['knotsU'], j['knotsV'], j['weights']),
  GridSurface: (j) => new GridSurface(grid(j['points'])),
  BilinearSurface: (j) => {
    const [a, b, c, d] = (j['corners'] as number[][]).map(point) as [Vector3, Vector3, Vector3, Vector3];
    return new BilinearSurface(a, b, c, d);
  },
  HermiteSurface: (j) =>
    new HermiteSurface({
      points: grid(j['points']),
      tangentsU: grid(j['tangentsU']),
      tangentsV: grid(j['tangentsV']),
      twists: grid(j['twists']),
      knotsU: j['knotsU'],
      knotsV: j['knotsV'],
    }),
  ExtrusionSurface: (j) => new ExtrusionSurface(curveFromJSON(j['curve']), point(j['direction']), j['v'][0], j['v'][1]),
  RevolutionSurface: (j) => new RevolutionSurface(curveFromJSON(j['curve']), point(j['axisOrigin']), point(j['axisDirection']), j['u'][0], j['u'][1]),
  RuledSurface: (j) => new RuledSurface(curveFromJSON(j['curveA']), curveFromJSON(j['curveB'])),
  LoftSurface: (j) => new LoftSurface(j['curves'].map(curveFromJSON), { parameters: j['parameters'], closed: j['closed'] }),
  CoonsSurface: (j) => new CoonsSurface(curveFromJSON(j['bottom']), curveFromJSON(j['top']), curveFromJSON(j['left']), curveFromJSON(j['right'])),
  GordonSurface: (j) =>
    new GordonSurface(j['curvesU'].map(curveFromJSON), j['curvesV'].map(curveFromJSON), { parametersU: j['parametersU'], parametersV: j['parametersV'] }),
  SweepSurface: (j) => new SweepSurface(curveFromJSON(j['profile']), curveFromJSON(j['path']), { mode: j['mode'] }),
  PipeSurface: (j) => new PipeSurface(curveFromJSON(j['path']), j['radius']),
  CompositeSurface: (j) => new CompositeSurface(j['patches'].map((column: SurfaceJSON[]) => column.map(surfaceFromJSON))),
  TrimmedSurface: (j) => new TrimmedSurface(surfaceFromJSON(j['surface']), j['u'][0], j['u'][1], j['v'][0], j['v'][1]),
  ReversedSurface: (j) => new ReversedSurface(surfaceFromJSON(j['surface'])),
  TransformedSurface: (j) => new TransformedSurface(surfaceFromJSON(j['surface']), new Transform3(j['linear'], point(j['translation']))),
  OffsetSurface: (j) => new OffsetSurface(surfaceFromJSON(j['surface']), j['distance']),
};

/** Rebuilds a surface from the output of `toJSON`, also after a pass through `JSON.stringify`. */
export function surfaceFromJSON(json: SurfaceJSON): ISurface {
  const reader = readers[json.type];
  if (!reader) {
    throw new TypeError(`Unknown surface type: ${json.type}`);
  }
  return reader(json);
}
