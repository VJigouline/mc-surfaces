import { BezierCurve, BSplineCurve, CircularArc, EllipticalArc, HyperbolicArc, Line, ParametricCurve, Polyline } from 'mc-curves';
import { describe, expect, it } from 'vitest';
import {
  BezierSurface,
  BilinearSurface,
  BSplineSurface,
  CompositeSurface,
  ConeSurface,
  CoonsSurface,
  curveFromJSON,
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
  IsoCurve,
  LoftSurface,
  NurbsSurface,
  OffsetSurface,
  ParametricSurface,
  PipeSurface,
  PlaneSurface,
  ReversedSurface,
  RevolutionSurface,
  RotationMinimizingFrames,
  RuledSurface,
  SphereSurface,
  surfaceFromJSON,
  SweepSurface,
  TorusSurface,
  Transform3,
  TransformedSurface,
  TrimmedSurface,
  Vector3,
} from '../src/index.js';
import { CardinalSplines } from '../src/math/cardinal.js';
import { tilted } from './catalog.js';
import { expectPoint, v } from './helpers.js';

const PI = Math.PI;

describe('PlaneSurface', () => {
  const plane = new PlaneSurface(tilted, { u: [-1, 2], v: [0, 4] });

  it('has its closed-form measures', () => {
    expect(plane.area).toBeCloseTo(12, 12);
    expect(plane.isPlanar()).toBe(true);
    expectPoint(plane.normalAt(0, 0), tilted.planeNormal);
    expectPoint(plane.planeNormal(), tilted.planeNormal);
    expect(plane.gaussianCurvatureAt(0, 1)).toBeCloseTo(0, 12);
    expect(plane.enclosesVolume).toBe(false);
    expect(plane.enclosedVolume()).toBeUndefined();
  });

  it('finds the closest point by projection, clamped to the rectangle', () => {
    const inside = plane.closestPoint(plane.pointAt(0.5, 1).addScaled(plane.normalAt(0, 0), 3));
    expect(inside.distance).toBeCloseTo(3, 12);
    expect(inside.u).toBeCloseTo(0.5, 12);
    const corner = plane.closestPoint(plane.pointAt(5, 9));
    expectPoint(corner.point, plane.pointAt(2, 4));
  });

  it('builds a centred patch, and scales keeping its parameters', () => {
    const patch = PlaneSurface.centred(v(1, 1, 1), Vector3.Z, 4, 2);
    expect(patch.area).toBeCloseTo(8, 12);
    expectPoint(patch.centroid(), v(1, 1, 1));
    const big = patch.transformed(Transform3.scale(3)) as PlaneSurface;
    expect(big).toBeInstanceOf(PlaneSurface);
    expect(big.area).toBeCloseTo(72, 10);
    expect(big.domainU.end).toBe(2);
  });

  it('has exact offset, iso-curves and NURBS form', () => {
    const moved = plane.offset(2);
    expect(moved).toBeInstanceOf(PlaneSurface);
    expectPoint(moved.pointAt(1, 1), plane.pointAt(1, 1).addScaled(tilted.planeNormal, 2));
    expect(plane.isoCurveU(1)).toBeInstanceOf(Line);
    expect(plane.isoCurveV(1)).toBeInstanceOf(Line);
    const nurbs = plane.toNurbs();
    expect(nurbs.degreeU).toBe(1);
    expectPoint(nurbs.pointAt(0.3, 2.2), plane.pointAt(0.3, 2.2));
  });

  it('rejects an empty domain', () => {
    expect(() => new PlaneSurface(tilted, { u: [1, 1], v: [0, 1] })).toThrow(RangeError);
    expect(() => new PlaneSurface(tilted, { u: [0, 1], v: [0, NaN] })).toThrow(RangeError);
    expect(() => new PlaneSurface(tilted, { u: [0, 1], v: [0, 1] }, 0)).toThrow(RangeError);
  });
});

describe('CylinderSurface', () => {
  const r = 2;
  const cylinder = new CylinderSurface(Frame3.WORLD, r, { v: [0, 5] });

  it('has its closed-form measures', () => {
    expect(cylinder.area).toBeCloseTo(2 * PI * r * 5, 10);
    expect(cylinder.isPeriodicU).toBe(true);
    expect(cylinder.isClosedU).toBe(true);
    expect(cylinder.isClosedV).toBe(false);
    expect(cylinder.enclosesVolume).toBe(false);
    expect(cylinder.isPlanar()).toBe(false);
  });

  it('agrees with the general area integrator', () => {
    const trimmed = new TrimmedSurface(cylinder, 0.5, 2, 1, 3);
    expect(trimmed.area).toBeCloseTo(r * 1.5 * 2, 8);
  });

  it('has an outward normal, and principal curvatures 0 and -1/r', () => {
    expectPoint(cylinder.normalAt(0, 1), Vector3.X);
    expectPoint(cylinder.normalAt(PI / 2, 1), Vector3.Y);
    const p = cylinder.principalCurvaturesAt(1, 2);
    expect(p.k1).toBeCloseTo(0, 10);
    expect(p.k2).toBeCloseTo(-1 / r, 10);
    expect(Math.abs(p.direction1.dot(Vector3.Z))).toBeCloseTo(1, 9);
    expect(cylinder.gaussianCurvatureAt(1, 2)).toBeCloseTo(0, 10);
    expect(cylinder.meanCurvatureAt(1, 2)).toBeCloseTo(-1 / (2 * r), 10);
  });

  it('finds closest points in closed form, also on a patch', () => {
    const c = cylinder.closestPoint(v(5, 0, 2));
    expect(c.distance).toBeCloseTo(3, 12);
    expectPoint(c.point, v(2, 0, 2));
    const patch = new CylinderSurface(Frame3.WORLD, r, { u: [0, PI / 2], v: [0, 1] });
    expectPoint(patch.closestPoint(v(-1, -0.1, 0.5)).point, patch.pointAt(0, 0.5).sub(v(0, 0, 0)), 3);
    const near = patch.closestPoint(v(1, -3, 7));
    expectPoint(near.point, patch.pointAt(0, 1));
  });

  it('offsets exactly, and falls back when the radius would vanish', () => {
    const fat = cylinder.offset(1) as CylinderSurface;
    expect(fat).toBeInstanceOf(CylinderSurface);
    expect(fat.radius).toBe(3);
    expect(cylinder.offset(-5)).toBeInstanceOf(OffsetSurface);
  });

  it('turns inside out in a mirrored frame, and offsets accordingly', () => {
    const mirrored = cylinder.transformed(Transform3.mirror(Vector3.Z)) as CylinderSurface;
    expect(mirrored).toBeInstanceOf(CylinderSurface);
    expectPoint(mirrored.normalAt(0, 1), Vector3.X.negate());
    expect(mirrored.offset(1).distanceTo(v(1, 0, -1))).toBeCloseTo(0, 9);
  });

  it('has an exact NURBS form', () => {
    const nurbs = cylinder.toNurbs();
    for (const p of nurbs.sample(11, 5).flat()) {
      expect(Math.hypot(p.x, p.y)).toBeCloseTo(r, 12);
    }
    expect(nurbs.area / cylinder.area).toBeCloseTo(1, 7);
  });

  it('turns under a non-uniform scale into a general surface', () => {
    const stretched = cylinder.transformed(Transform3.scale(1, 2, 1));
    expect(stretched).toBeInstanceOf(TransformedSurface);
    expectPoint(stretched.pointAt(PI / 2, 1), v(0, 4, 1));
    expect(stretched.toNurbs().distanceTo(v(0, 4, 2))).toBeLessThan(1e-9);
  });
});

describe('ConeSurface', () => {
  const angle = PI / 6;
  const cone = new ConeSurface(Frame3.WORLD, 0, angle, { v: [0, 3] });

  it('has the area of a cone', () => {
    const radius = 3 * Math.tan(angle);
    const slant = 3 / Math.cos(angle);
    expect(cone.area).toBeCloseTo(PI * radius * slant, 10);
    expect(new TrimmedSurface(cone, 0, 2 * PI, 0, 3).area / cone.area).toBeCloseTo(1, 8);
  });

  it('has a limit normal at the apex', () => {
    const n = cone.normalAt(0, 0);
    expect(n.length).toBeCloseTo(1, 9);
    expectPoint(n, cone.normalAt(0, 1e-3), 1e-6);
  });

  it('has curvature zero and -cos(a)/r', () => {
    const p = cone.principalCurvaturesAt(1, 2);
    const r = 2 * Math.tan(angle);
    expect(p.k1).toBeCloseTo(0, 9);
    expect(p.k2).toBeCloseTo(-Math.cos(angle) / r, 9);
  });

  it('rejects a radius going negative or a right angle', () => {
    expect(() => new ConeSurface(Frame3.WORLD, 1, -PI / 4, { v: [0, 3] })).toThrow(RangeError);
    expect(() => new ConeSurface(Frame3.WORLD, 1, PI / 2)).toThrow(RangeError);
    expect(() => new ConeSurface(Frame3.WORLD, 0, 0)).toThrow(RangeError);
  });

  it('offsets to a cone with the same angle', () => {
    const frustum = new ConeSurface(tilted, 1, angle, { v: [1, 3] });
    const shell = frustum.offset(0.5) as ConeSurface;
    expect(shell).toBeInstanceOf(ConeSurface);
    expect(shell.semiAngle).toBe(angle);
    const target = frustum.pointAt(1, 2).addScaled(frustum.normalAt(1, 2), 0.5);
    expect(shell.distanceTo(target)).toBeLessThan(1e-7);
    expect(new ConeSurface(Frame3.WORLD, 0.1, angle).offset(-1)).toBeInstanceOf(OffsetSurface);
    // A wrong shift of v stays on the same infinite cone, so check the patch's own edges.
    const lower = frustum.pointAt(1, 1).addScaled(frustum.normalAt(1, 1), 0.5);
    const upper = frustum.pointAt(1, 3).addScaled(frustum.normalAt(1, 3), 0.5);
    expectPoint(shell.pointAt(1, shell.domainV.start), lower, 1e-9);
    expectPoint(shell.pointAt(1, shell.domainV.end), upper, 1e-9);
  });
});

describe('SphereSurface', () => {
  const r = 2.5;
  const sphere = SphereSurface.centredAt(v(1, 2, 3), r);

  it('has the area and volume of a sphere', () => {
    expect(sphere.area).toBeCloseTo(4 * PI * r * r, 10);
    expect(sphere.enclosesVolume).toBe(true);
    expect(sphere.enclosedVolume()).toBeCloseTo((4 / 3) * PI * r ** 3, 7);
    expectPoint(sphere.centroid(), v(1, 2, 3), 1e-7);
  });

  it('has the area of a zone', () => {
    const zone = new SphereSurface(Frame3.WORLD, r, { v: [0, PI / 2] });
    expect(zone.area).toBeCloseTo(2 * PI * r * r, 10);
    expect(zone.enclosesVolume).toBe(false);
  });

  it('has an outward normal everywhere, also at the poles', () => {
    for (const [a, b] of [[0.3, 0.2], [2, -1], [5, 1.4]] as const) {
      expectPoint(sphere.normalAt(a, b), sphere.pointAt(a, b).sub(v(1, 2, 3)).normalized());
    }
    expectPoint(sphere.normalAt(1, PI / 2), Vector3.Z);
    expectPoint(sphere.normalAt(1, -PI / 2), Vector3.Z.negate());
  });

  it('has curvature -1/r everywhere, also at the poles', () => {
    for (const [a, b] of [[0.3, 0.2], [2, -1], [1, PI / 2], [1, -PI / 2]] as const) {
      const p = sphere.principalCurvaturesAt(a, b);
      expect(p.k1).toBeCloseTo(-1 / r, 6);
      expect(p.k2).toBeCloseTo(-1 / r, 6);
      expect(sphere.gaussianCurvatureAt(a, b)).toBeCloseTo(1 / (r * r), 5);
    }
  });

  it('gives usable tangents at a pole', () => {
    const f = sphere.frameAt(0.7, PI / 2);
    expect(f.tangentU.length).toBeCloseTo(1, 9);
    expect(f.tangentV.length).toBeCloseTo(1, 9);
    expect(f.tangentU.dot(f.normal)).toBeCloseTo(0, 3);
  });

  it('wraps longitude and clamps latitude', () => {
    expectPoint(sphere.pointAt(0.5 + 2 * PI, 0.2), sphere.pointAt(0.5, 0.2));
    expectPoint(sphere.pointAt(0.5, 3), sphere.pointAt(0.5, PI / 2));
  });

  it('finds closest points directly, or on the edge of a patch', () => {
    const c = sphere.closestPoint(v(1, 2, 13));
    expect(c.distance).toBeCloseTo(10 - r, 9);
    expectPoint(c.point, v(1, 2, 3 + r));
    const patch = new SphereSurface(Frame3.WORLD, r, { u: [0, 1], v: [0, 1] });
    const edge = patch.closestPoint(v(-5, -5, -5));
    expect(edge.u === 0 || edge.u === 1 || Math.abs(edge.v) < 1e-9).toBe(true);
  });

  it('has an exact NURBS form', () => {
    const nurbs = sphere.toNurbs();
    expect(nurbs.isRational).toBe(true);
    for (const p of nurbs.sample(21, 21).flat()) {
      expect(p.distanceTo(v(1, 2, 3))).toBeCloseTo(r, 11);
    }
    expect(nurbs.area / sphere.area).toBeCloseTo(1, 6);
    expect(nurbs.enclosedVolume()! / sphere.enclosedVolume()!).toBeCloseTo(1, 5);
  });

  it('offsets by changing the radius, keeps the class when scaled', () => {
    expect((sphere.offset(1) as SphereSurface).radius).toBe(3.5);
    expect(sphere.offset(-r - 1)).toBeInstanceOf(OffsetSurface);
    expect((sphere.transformed(Transform3.scale(2)) as SphereSurface).radius).toBe(2 * r);
  });

  it('rejects a latitude beyond the poles and a bad radius', () => {
    expect(() => new SphereSurface(Frame3.WORLD, 1, { v: [0, 2] })).toThrow(RangeError);
    expect(() => new SphereSurface(Frame3.WORLD, -1)).toThrow(RangeError);
  });
});

describe('TorusSurface', () => {
  const [R, r] = [3, 1];
  const torus = new TorusSurface(tilted, R, r);

  it('has the area and volume of a torus', () => {
    expect(torus.area).toBeCloseTo(4 * PI * PI * R * r, 9);
    expect(torus.enclosesVolume).toBe(true);
    expect(torus.enclosedVolume()).toBeCloseTo(2 * PI * PI * R * r * r, 7);
    expectPoint(torus.centroid(), tilted.origin, 1e-7);
  });

  it('has its known curvature on the outside and inside', () => {
    expect(torus.gaussianCurvatureAt(1, 0)).toBeCloseTo(1 / (r * (R + r)), 9);
    expect(torus.gaussianCurvatureAt(1, PI)).toBeCloseTo(-1 / (r * (R - r)), 9);
    expect(torus.gaussianCurvatureAt(1, PI / 2)).toBeCloseTo(0, 9);
  });

  it('is periodic both ways', () => {
    expect(torus.isPeriodicU && torus.isPeriodicV).toBe(true);
    expectPoint(torus.pointAt(1 + 2 * PI, 1 - 2 * PI), torus.pointAt(1, 1), 1e-12);
  });

  it('has an exact NURBS form and offset', () => {
    const nurbs = torus.toNurbs();
    for (const p of nurbs.sample(13, 13).flat()) {
      const local = tilted.toLocal(p);
      expect(Math.hypot(Math.hypot(local.x, local.y) - R, local.z)).toBeCloseTo(r, 10);
    }
    expect((torus.offset(0.5) as TorusSurface).minorRadius).toBe(1.5);
    expect(torus.offset(-2)).toBeInstanceOf(OffsetSurface);
  });
});

describe('quadrics', () => {
  it('has an ellipsoid that reduces to a sphere, with the right volume', () => {
    const ball = new EllipsoidSurface(Frame3.WORLD, 2, 2, 2);
    expect(ball.area).toBeCloseTo(16 * PI, 8);
    const e = new EllipsoidSurface(tilted, 3, 2, 1);
    expect(e.enclosedVolume()).toBeCloseTo((4 / 3) * PI * 6, 7);
    expect(() => new EllipsoidSurface(tilted, 0, 1, 1)).toThrow(RangeError);
    expect(() => new EllipsoidSurface(tilted, 1, 1, 1, { v: [0, 2] })).toThrow(RangeError);
  });

  it('has an exact NURBS ellipsoid', () => {
    const nurbs = new EllipsoidSurface(Frame3.WORLD, 3, 2, 1).toNurbs();
    for (const p of nurbs.sample(17, 17).flat()) {
      expect((p.x / 3) ** 2 + (p.y / 2) ** 2 + p.z ** 2).toBeCloseTo(1, 11);
    }
  });

  it('has a paraboloid meeting its equation, with a limit normal at the vertex', () => {
    const para = new EllipticParaboloidSurface(Frame3.WORLD, 2, 1, 0.5);
    const p = para.pointAt(0.7, 0.8);
    expect(p.z).toBeCloseTo(0.5 * ((p.x / 2) ** 2 + p.y ** 2), 12);
    expectPoint(para.normalAt(1, 0), Vector3.Z.negate(), 1e-9);
    expect(() => new EllipticParaboloidSurface(Frame3.WORLD, 1, 1, 1, { v: [-1, 1] })).toThrow(RangeError);
    for (const q of para.toNurbs().sample(13, 13).flat()) {
      expect(q.z).toBeCloseTo(0.5 * ((q.x / 2) ** 2 + q.y ** 2), 11);
    }
  });

  it('has a saddle with negative curvature and a matching Bezier form', () => {
    const saddle = new HyperbolicParaboloidSurface(Frame3.WORLD, 1, 1, { u: [-1, 1], v: [-1, 1] });
    expect(saddle.pointAt(0.5, 0.25).z).toBeCloseTo(0.25 - 0.0625, 12);
    expect(saddle.gaussianCurvatureAt(0, 0)).toBeCloseTo(-4, 9);
    const nurbs = saddle.toNurbs();
    expect(nurbs).toBeInstanceOf(BezierSurface);
    expectPoint(nurbs.pointAt(0.75, 0.25), saddle.pointAt(0.5, -0.5), 1e-12);
    const big = saddle.transformed(Transform3.scale(2));
    expectPoint(big.pointAt(0.5, 0.25), saddle.pointAt(0.5, 0.25).scale(2));
  });

  it('has hyperboloids meeting their equations, in both forms', () => {
    const one = new HyperboloidSurface(Frame3.WORLD, 2, 1, 3);
    const two = new HyperboloidSurface(Frame3.WORLD, 2, 1, 3, { v: [0, 1.5] }, 2);
    const eq = (p: Vector3) => (p.x / 2) ** 2 + p.y ** 2 - (p.z / 3) ** 2;
    expect(eq(one.pointAt(1, 0.5))).toBeCloseTo(1, 12);
    expect(eq(two.pointAt(1, 0.5))).toBeCloseTo(-1, 12);
    for (const q of one.toNurbs().sample(9, 9).flat()) {
      expect(eq(q)).toBeCloseTo(1, 10);
    }
    for (const q of two.toNurbs().sample(9, 9).flat()) {
      expect(eq(q)).toBeCloseTo(-1, 10);
    }
    expect(() => new HyperboloidSurface(Frame3.WORLD, 1, 1, 1, { v: [-1, 1] }, 2)).toThrow(RangeError);
    expect(() => new HyperbolicArc(Frame3.WORLD, 1, 1, 0, 1).length).not.toThrow();
  });

  it('has a helicoid made of straight lines, and is a minimal surface', () => {
    const screw = new HelicoidSurface(Frame3.WORLD, 2, { u: [0, 4 * PI], v: [0.5, 2] });
    expectPoint(screw.pointAt(PI / 2, 1), v(0, 1, 0.5));
    expect(screw.meanCurvatureAt(1, 1)).toBeCloseTo(0, 9);
    expect(screw.isPeriodicU).toBe(false);
    expect(() => new HelicoidSurface(Frame3.WORLD, 0, { u: [0, 1], v: [0, 1] })).toThrow(RangeError);
    expect(screw.degreeV).toBe(1);
  });
});

describe('NurbsSurface', () => {
  const net = Array.from({ length: 5 }, (_, i) => Array.from({ length: 4 }, (_, j) => v(i, j, Math.sin(i) * Math.cos(j))));
  const spline = BSplineSurface.uniform(3, 2, net);

  it('interpolates its corners when clamped', () => {
    expectPoint(spline.pointAt(0, 0), net[0]![0]!);
    expectPoint(spline.pointAt(1, 1), net[4]![3]!);
    expect(spline.type).toBe('BSplineSurface');
    expect(spline.degreeU).toBe(3);
    expect(spline.degreeV).toBe(2);
    expect(spline.breakpointsU()).toEqual([0, 0.5, 1]);
  });

  it('forms a partition of unity, so a flat net stays flat', () => {
    const flat = BSplineSurface.uniform(2, 2, net.map((row) => row.map((p) => v(p.x, p.y, 7))));
    for (const q of flat.sample(9, 9).flat()) {
      expect(q.z).toBeCloseTo(7, 12);
    }
    expect(flat.isPlanar()).toBe(true);
    expect(Math.abs(flat.planeNormal()!.z)).toBeCloseTo(1, 12);
  });

  it('keeps its shape through knot insertion', () => {
    const more = spline.insertKnotU(0.3, 2).insertKnotV(0.6);
    expect(more.controlPoints.length).toBe(net.length + 2);
    expect(more.controlPoints[0]!.length).toBe(net[0]!.length + 1);
    for (const [a, b] of [[0.1, 0.2], [0.3, 0.6], [0.77, 0.9]] as const) {
      expectPoint(more.pointAt(a, b), spline.pointAt(a, b), 1e-12);
    }
  });

  it('keeps a rational shape through trimming, including equal weights on a line', () => {
    const rational = new NurbsSurface(2, 2, net.slice(0, 3).map((row) => row.slice(0, 3)), [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1], [
      [2, 2, 2],
      [1, 3, 1],
      [2, 2, 2],
    ]);
    const piece = rational.trim(0.2, 0.7, 0.1, 0.8);
    expect(piece).toBeInstanceOf(NurbsSurface);
    expectPoint(piece.pointAt(0.5, 0.5), rational.pointAt(0.5, 0.5), 1e-11);
    expectPoint(piece.pointAt(0.2, 0.1), rational.pointAt(0.2, 0.1), 1e-11);
  });

  it('keeps a column of equal weights through trimming', () => {
    const grid = net.slice(0, 3).map((row) => row.slice(0, 3));
    const rational = new NurbsSurface(2, 2, grid, [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1], [
      [2, 1, 1],
      [2, 3, 1],
      [2, 1, 4],
    ]);
    const piece = rational.trim(0.2, 0.7, 0.1, 0.8);
    for (const [a, b] of [[0.25, 0.15], [0.5, 0.5], [0.65, 0.75]] as const) {
      expectPoint(piece.pointAt(a, b), rational.pointAt(a, b), 1e-11);
    }
  });

  it('has exact iso-curves, also when rational', () => {
    const rational = new NurbsSurface(2, 2, net.slice(0, 3).map((row) => row.slice(0, 3)), [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1], [
      [1, 2, 1],
      [3, 1, 2],
      [1, 1, 4],
    ]);
    for (const s of [rational, spline]) {
      expectPoint(s.isoCurveU(0.3).pointAt(0.8), s.pointAt(0.3, 0.8), 1e-12);
      expectPoint(s.isoCurveV(0.8).pointAt(0.3), s.pointAt(0.3, 0.8), 1e-12);
    }
  });

  it('computes rational derivatives through the quotient rule', () => {
    const quarter = new SphereSurface(Frame3.WORLD, 1, { u: [0, PI / 2], v: [0, PI / 2] }).toNurbs();
    const n = quarter.normalAt(quarter.domainU.start + 0.4, 0.3);
    const p = quarter.pointAt(quarter.domainU.start + 0.4, 0.3);
    expectPoint(n, p.normalized(), 1e-9);
    expect(quarter.gaussianCurvatureAt(0.4, 0.3)).toBeCloseTo(1, 8);
  });

  it('interpolates a grid', () => {
    const smooth = BSplineSurface.interpolate(net);
    const us = smooth.domainU;
    expect(us.start).toBe(0);
    for (const row of net) {
      for (const p of row) {
        expect(smooth.distanceTo(p)).toBeLessThan(1e-8);
      }
    }
    expect(() => BSplineSurface.interpolate([[v(0, 0)]])).toThrow(RangeError);
    expect(() => BSplineSurface.interpolate(net.map((row, i) => (i === 2 ? row.slice(1) : row)))).toThrow(RangeError);
  });

  it('rejects invalid definitions', () => {
    const knots = [0, 0, 1, 1];
    const grid = [[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1)]];
    expect(() => new NurbsSurface(0, 1, grid, knots, knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, grid, [0, 1], knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, grid, [1, 0, 1, 1], knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, grid, [0, 0, 0, 0], knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(2, 1, grid, [0, 0, 0, 1, 1], knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, [[v(0, 0), v(0, 1)], [v(1, 0)]], knots, knots)).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, grid, knots, knots, [[1, 1]])).toThrow(RangeError);
    expect(() => new NurbsSurface(1, 1, grid, knots, knots, [[1, 0], [1, 1]])).toThrow(RangeError);
    expect(() => BSplineSurface.uniform(3, 3, grid)).toThrow(RangeError);
  });

  it('treats equal weights as plain', () => {
    const grid = [[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1)]];
    const same = new NurbsSurface(1, 1, grid, [0, 0, 1, 1], [0, 0, 1, 1], [[2, 2], [2, 2]]);
    expect(same.isRational).toBe(false);
    expect(same.type).toBe('BSplineSurface');
  });
});

describe('BezierSurface', () => {
  it('matches the Bernstein form', () => {
    const net = [[v(0, 0), v(0, 1, 2)], [v(1, 0, 1), v(1, 1)], [v(2, 0), v(2, 1, 3)]];
    const s = new BezierSurface(net);
    const [a, b] = [0.3, 0.6];
    const bu = [(1 - a) ** 2, 2 * a * (1 - a), a * a];
    const bv = [1 - b, b];
    let sum = Vector3.ZERO;
    net.forEach((row, i) => row.forEach((p, j) => {
      sum = sum.addScaled(p, bu[i]! * bv[j]!);
    }));
    expectPoint(s.pointAt(a, b), sum, 1e-12);
    expect(s.degreeU).toBe(2);
    expect(s.degreeV).toBe(1);
    expect(s.transformed(Transform3.translate(v(1, 1, 1)))).toBeInstanceOf(BezierSurface);
    expect(() => new BezierSurface([[v(0, 0)]])).toThrow(RangeError);
  });

  it('reproduces a sphere patch when rational', () => {
    const quarterArc = new CircularArc(Frame3.WORLD, 1, 0, PI / 2).toNurbs();
    expect(quarterArc.isRational).toBe(true);
    const w = Math.SQRT1_2;
    const octant = new BezierSurface(
      [
        [v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)],
        [v(1, 1, 0), v(1, 1, 1), v(0, 0, 1)],
        [v(0, 1, 0), v(0, 1, 1), v(0, 0, 1)],
      ],
      [
        [1, w, 1],
        [w, w * w, w],
        [1, w, 1],
      ],
    );
    expect(octant.type).toBe('RationalBezierSurface');
    for (const p of octant.sample(9, 9).flat()) {
      expect(p.length).toBeCloseTo(1, 12);
    }
    expect(octant.area).toBeCloseTo(PI / 2, 7);
  });
});

describe('patches', () => {
  it('has a bilinear patch through its corners', () => {
    const b = new BilinearSurface(v(0, 0), v(2, 0), v(0, 2), v(2, 2, 2));
    expectPoint(b.pointAt(1, 0), v(2, 0));
    expectPoint(b.pointAt(0, 1), v(0, 2));
    expectPoint(b.pointAt(0.5, 0.5), v(1, 1, 0.5));
    expect(b.type).toBe('BilinearSurface');
    expect(b.transformed(Transform3.scale(2))).toBeInstanceOf(BilinearSurface);
    const flat = new BilinearSurface(v(0, 0), v(2, 0), v(0, 3), v(2, 3));
    expect(flat.area).toBeCloseTo(6, 12);
  });

  it('has a grid addressed by index, right-handed at joints', () => {
    const g = new GridSurface([[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1)], [v(3, 0), v(3, 1)]]);
    expectPoint(g.pointAt(1.5, 0.5), v(2, 0.5));
    expectPoint(g.derivativeAt(1, 0.5, 1, 0), v(2, 0));
    expectPoint(g.derivativeAt(0.99, 0.5, 1, 0), v(1, 0));
    expect(g.derivativeAt(0.5, 0.5, 2, 0).length).toBe(0);
    expect(g.breakpointsU()).toEqual([0, 1, 2]);
    expect(g.area).toBeCloseTo(3, 12);
    expectPoint(g.toNurbs().pointAt(1.5, 0.5), v(2, 0.5), 1e-12);
    expect(() => new GridSurface([[v(0, 0)]])).toThrow(RangeError);
    expect(() => new GridSurface([[v(0, 0), v(1, 0)], [v(0, 1)]])).toThrow(RangeError);
  });

  it('has a Hermite surface that matches its data', () => {
    const data = {
      points: [[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1, 1)]],
      tangentsU: [[v(1, 0), v(1, 0, 1)], [v(1, 0, 2), v(1, 0)]],
      tangentsV: [[v(0, 1), v(0, 1, 3)], [v(0, 1), v(0, 1, -1)]],
    };
    const h = new HermiteSurface(data);
    expectPoint(h.pointAt(1, 1), v(1, 1, 1));
    expectPoint(h.derivativeAt(1, 0, 1, 0), v(1, 0, 2));
    expectPoint(h.derivativeAt(0, 1, 0, 1), v(0, 1, 3));
    expectPoint(h.derivativeAt(0, 0, 1, 1), Vector3.ZERO);
    expect(h.degreeU).toBe(3);
    expect(() => new HermiteSurface({ ...data, tangentsU: [] })).toThrow(RangeError);
    expect(() => new HermiteSurface({ ...data, knotsU: [0, 0] })).toThrow(RangeError);
    expect(() => new HermiteSurface({ ...data, knotsU: [0] })).toThrow(RangeError);
    expect(h.derivativeAt(0.5, 0.5, 4, 0).length).toBe(0);
  });

  it('takes the twist given at each corner', () => {
    const data = {
      points: [[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1)]],
      tangentsU: [[v(1, 0), v(1, 0)], [v(1, 0), v(1, 0)]],
      tangentsV: [[v(0, 1), v(0, 1)], [v(0, 1), v(0, 1)]],
      twists: [[v(0, 0, 1), v(0, 0, -2)], [v(0, 0, 3), v(0, 0, 4)]],
    };
    const h = new HermiteSurface(data);
    expectPoint(h.derivativeAt(0, 0, 1, 1), v(0, 0, 1), 1e-12);
    expectPoint(h.derivativeAt(0, 0.999999, 1, 1), v(0, 0, -2), 1e-4);
    expectPoint(h.derivativeAt(0.999999, 0, 1, 1), v(0, 0, 3), 1e-4);
    expect(Math.abs(h.pointAt(0.5, 0.5).z)).toBeGreaterThan(0.01);
  });

  it('converts a Hermite surface to NURBS exactly', () => {
    const data = {
      points: [[v(0, 0), v(0, 1)], [v(1, 0), v(1, 1, 1)]],
      tangentsU: [[v(1, 0), v(1, 0, 1)], [v(1, 0, 2), v(1, 0)]],
      tangentsV: [[v(0, 1), v(0, 1, 3)], [v(0, 1), v(0, 1, -1)]],
      twists: [[v(0, 0, 1), v(0, 0, -1)], [v(0, 0, 0.5), v(0, 0, 2)]],
    };
    const h = new HermiteSurface(data);
    const n = h.toNurbs(1e-3);
    for (const [a, b] of [[0.1, 0.2], [0.5, 0.5], [0.77, 0.31]] as const) {
      expectPoint(n.pointAt(a, b), h.pointAt(a, b), 1e-11);
    }
  });
});

describe('tessellation', () => {
  it('refines for twist alone, where every edge is straight', () => {
    const twisted = new BilinearSurface(v(0, 0), v(1, 0), v(0, 1), v(1, 1, 8));
    const chordTolerance = 1e-3;
    const mesh = twisted.tessellate({ chordTolerance });
    for (let t = 0; t < mesh.triangles.length; t += 3) {
      const [a, b] = [mesh.parameters[mesh.triangles[t]!]!, mesh.parameters[mesh.triangles[t + 1]!]!];
      const c = mesh.parameters[mesh.triangles[t + 2]!]!;
      const flat = mesh.points[mesh.triangles[t]!]!.lerp(mesh.points[mesh.triangles[t + 2]!]!, 0.5);
      const onSurface = twisted.pointAt((a.u + c.u) / 2, (a.v + c.v) / 2);
      expect(flat.distanceTo(onSurface)).toBeLessThan(1.5 * chordTolerance);
      expect(b).toBeDefined();
    }
  });

  it('honours the edge length and point limits', () => {
    const plate = new PlaneSurface(Frame3.WORLD, { u: [0, 10], v: [0, 10] });
    const fine = plate.tessellate({ maxEdgeLength: 0.5 });
    for (let t = 0; t < fine.triangles.length; t += 3) {
      const [a, b] = [fine.points[fine.triangles[t]!]!, fine.points[fine.triangles[t + 1]!]!];
      expect(a.distanceTo(b)).toBeLessThanOrEqual(0.5 + 1e-9);
    }
    const sphere = SphereSurface.centredAt(Vector3.ZERO, 1);
    expect(sphere.tessellate({ chordTolerance: 1e-6, maxPoints: 2000 }).points.length).toBeLessThanOrEqual(2000);
  });
});

describe('surfaces from curves', () => {
  it('has an extrusion of a circle equal to a cylinder', () => {
    const e = new ExtrusionSurface(CircularArc.circle(v(0, 0), 2), v(0, 0, 5));
    expect(e.area).toBeCloseTo(2 * PI * 2 * 5, 9);
    expect(e.isPeriodicU).toBe(true);
    const nurbs = e.toNurbs();
    expect(nurbs.isRational).toBe(true);
    for (const p of nurbs.sample(9, 3).flat()) {
      expect(Math.hypot(p.x, p.y)).toBeCloseTo(2, 11);
    }
    expect(() => new ExtrusionSurface(new Line(v(0, 0), v(1, 0)), Vector3.ZERO)).toThrow(RangeError);
    expect(() => new ExtrusionSurface(new Line(v(0, 0), v(1, 0)), Vector3.Z, 1, 1)).toThrow(RangeError);
  });

  it('keeps an extrusion through any affine transform', () => {
    const e = new ExtrusionSurface(new Line(v(0, 0), v(1, 0)), Vector3.Z);
    const shear = new Transform3([1, 1, 0, 0, 1, 0, 0, 0, 2], v(0, 0, 0));
    const moved = e.transformed(shear);
    expect(moved).toBeInstanceOf(ExtrusionSurface);
    expectPoint(moved.pointAt(0.5, 0.5), shear.applyToPoint(e.pointAt(0.5, 0.5)));
    expect(moved.area).toBeGreaterThan(0);
  });

  it('has a revolution of a line equal to a cylinder, and of an arc equal to a sphere', () => {
    const tube = new RevolutionSurface(new Line(v(2, 0, 0), v(2, 0, 5)), Vector3.ZERO, Vector3.Z);
    const cylinder = new CylinderSurface(Frame3.WORLD, 2, { v: [0, 1] }, 5);
    expect(tube.equals(cylinder, 1e-12)).toBe(true);
    expect(tube.isPeriodicU).toBe(true);
    const halfFrame = Frame3.fromAxes(Vector3.ZERO, Vector3.X, Vector3.Z, Vector3.Y.negate());
    const ball = new RevolutionSurface(new CircularArc(halfFrame, 3, -PI / 2, PI / 2), Vector3.ZERO, Vector3.Z);
    expect(ball.area).toBeCloseTo(4 * PI * 9, 8);
    expect(ball.enclosesVolume).toBe(true);
    expect(ball.enclosedVolume()).toBeCloseTo((4 / 3) * PI * 27, 7);
    const nurbs = ball.toNurbs();
    for (const p of nurbs.sample(13, 13).flat()) {
      expect(p.length).toBeCloseTo(3, 11);
    }
  });

  it('has a part revolution, and handles profiles that touch the axis', () => {
    const cap = new RevolutionSurface(new Line(Vector3.ZERO, v(2, 0, 3)), Vector3.ZERO, Vector3.Z, 0, PI);
    expect(cap.isPeriodicU).toBe(false);
    expect(cap.area).toBeCloseTo(PI * 2 * Math.hypot(2, 3) * 0.5, 9);
    const nurbs = cap.toNurbs();
    expect(nurbs.distanceTo(cap.pointAt(1, 0.5))).toBeLessThan(1e-9);
  });

  it('keeps a revolution under rotation and scale, but not a mirror or shear', () => {
    const r = new RevolutionSurface(new Line(v(2, 0, 0), v(3, 0, 1)), Vector3.ZERO, Vector3.Z);
    expect(r.transformed(Transform3.rotate(Vector3.X, 1))).toBeInstanceOf(RevolutionSurface);
    expect(r.transformed(Transform3.mirror(Vector3.X))).toBeInstanceOf(TransformedSurface);
    expect(r.transformed(Transform3.scale(1, 2, 1))).toBeInstanceOf(TransformedSurface);
    expect(() => new RevolutionSurface(new Line(v(0, 0), v(1, 0)), Vector3.ZERO, Vector3.ZERO)).toThrow(RangeError);
    expect(() => new RevolutionSurface(new Line(v(0, 0), v(1, 0)), Vector3.ZERO, Vector3.Z, 1, 0)).toThrow(RangeError);
  });

  it('has a ruled surface of straight lines, whatever the curve domains', () => {
    const a = new CircularArc(Frame3.WORLD, 1, 0, PI);
    const b = new Line(v(-1, 3, 2), v(1, 3, 2));
    const r = new RuledSurface(a, b);
    expectPoint(r.pointAt(0.5, 0), a.pointAt(PI / 2));
    expectPoint(r.pointAt(0.5, 1), b.pointAt(0.5));
    expectPoint(r.derivativeAt(0.5, 0.3, 0, 1), b.pointAt(0.5).sub(a.pointAt(PI / 2)));
    expect(r.derivativeAt(0.5, 0.3, 0, 2).length).toBe(0);
    expect(r.degreeV).toBe(1);
    const flat = new RuledSurface(new Line(v(0, 0), v(2, 0)), new Line(v(0, 3), v(2, 3)));
    expect(flat.area).toBeCloseTo(6, 11);
  });

  it('has a loft that passes through its curves, and is C2 across them', () => {
    const curves = [0, 1, 3, 4].map((z, i) => new BezierCurve([v(0, 0, z), v(1, 1 + i * 0.3, z), v(2, 0, z)]));
    const loft = new LoftSurface(curves);
    expect(loft.parameters[0]).toBe(0);
    expect(loft.parameters[3]).toBe(1);
    loft.parameters.forEach((p, i) => {
      for (const s of [0, 0.4, 1]) {
        expectPoint(loft.pointAt(s, p), curves[i]!.pointAt(s), 1e-12);
      }
    });
    const joint = loft.parameters[1]!;
    const e = 1e-8;
    expectPoint(loft.derivativeAt(0.4, joint - e, 0, 1), loft.derivativeAt(0.4, joint + e, 0, 1), 1e-5);
    expectPoint(loft.derivativeAt(0.4, joint - e, 0, 2), loft.derivativeAt(0.4, joint + e, 0, 2), 1e-5);
    expect(loft.derivativeAt(0.4, 0.5, 0, 4).length).toBe(0);
  });

  it('has a loft of two curves equal to a ruled surface', () => {
    const [a, b] = [new Line(v(0, 0), v(2, 0)), new BezierCurve([v(0, 1, 1), v(1, 2, 1), v(2, 1, 1)])];
    expect(new LoftSurface([a, b]).equals(new RuledSurface(a, b), 1e-12)).toBe(true);
  });

  it('has a closed loft that returns smoothly to its first curve', () => {
    const curves = [0, 1, 2].map((k) => {
      const t = (k * 2 * PI) / 3;
      return new Line(v(Math.cos(t), Math.sin(t), 0), v(Math.cos(t), Math.sin(t), 1));
    });
    const loft = new LoftSurface(curves, { closed: true });
    expect(loft.isPeriodicV).toBe(true);
    expect(loft.isClosedV).toBe(true);
    expectPoint(loft.pointAt(0.3, 0), loft.pointAt(0.3, 1), 1e-12);
    expectPoint(loft.derivativeAt(0.3, 1e-9, 0, 1), loft.derivativeAt(0.3, 1 - 1e-9, 0, 1), 1e-6);
    expect(() => new LoftSurface(curves.slice(0, 2), { closed: true })).toThrow(RangeError);
  });

  it('rejects an invalid loft', () => {
    const a = new Line(v(0, 0), v(1, 0));
    const b = new Line(v(0, 1), v(1, 1));
    expect(() => new LoftSurface([a])).toThrow(RangeError);
    expect(() => new LoftSurface([a, a])).toThrow(RangeError);
    expect(() => new LoftSurface([a, b], { parameters: [0] })).toThrow(RangeError);
    expect(() => new LoftSurface([a, b], { parameters: [1, 0] })).toThrow(RangeError);
  });

  it('has a Coons patch that matches its four edges', () => {
    const bottom = new BezierCurve([v(0, 0, 0), v(1, -1, 1), v(2, 0, 0)]);
    const top = new BezierCurve([v(0, 2, 0), v(1, 3, 1), v(2, 2, 0)]);
    const left = new BezierCurve([v(0, 0, 0), v(-1, 1, 1), v(0, 2, 0)]);
    const right = new BezierCurve([v(2, 0, 0), v(3, 1, 1), v(2, 2, 0)]);
    const c = new CoonsSurface(bottom, top, left, right);
    for (const s of [0, 0.3, 0.7, 1]) {
      expectPoint(c.pointAt(s, 0), bottom.pointAt(s), 1e-12);
      expectPoint(c.pointAt(s, 1), top.pointAt(s), 1e-12);
      expectPoint(c.pointAt(0, s), left.pointAt(s), 1e-12);
      expectPoint(c.pointAt(1, s), right.pointAt(s), 1e-12);
    }
    expect(() => new CoonsSurface(new Line(v(9, 9), v(2, 0)), top, left, right)).toThrow(RangeError);
    const flat = new CoonsSurface(new Line(v(0, 0), v(2, 0)), new Line(v(0, 3), v(2, 3)), new Line(v(0, 0), v(0, 3)), new Line(v(2, 0), v(2, 3)));
    expect(flat.area).toBeCloseTo(6, 11);
    expect(flat.isPlanar()).toBe(true);
    expect(c.breakpointsU()).toEqual([0, 1]);
  });

  it('has a Gordon surface through all of its network', () => {
    // Each curve runs evenly in its coordinate, so they all cross at the same normalized parameters.
    const lift = (x: number, y: number) => 0.5 * Math.sin(x) * Math.cos(y);
    const uCurves = [0, 1, 2, 3].map((j) => new ParametricCurve((t) => v(2 * t, j, lift(2 * t, j)), 0, 1));
    const vCurves = [0, 1, 2].map((i) => new ParametricCurve((t) => v(i, 3 * t, lift(i, 3 * t)), 0, 1));
    const params = (n: number) => Array.from({ length: n }, (_, k) => k / (n - 1));
    const g = new GordonSurface(uCurves, vCurves, { parametersU: params(3), parametersV: params(4) });
    uCurves.forEach((curve, j) => {
      for (const s of [0.1, 0.5, 0.9]) {
        expectPoint(g.pointAt(s, j / 3), curve.pointAt(s), 1e-12);
      }
    });
    vCurves.forEach((curve, i) => {
      for (const s of [0.1, 0.5, 0.9]) {
        expectPoint(g.pointAt(i / 2, s), curve.pointAt(s), 1e-12);
      }
    });
  });

  it('finds the parameters of a network by itself', () => {
    const uCurves = [0, 1, 2].map((j) => new Line(v(0, j, 0), v(4, j, 0)));
    const vCurves = [0, 1, 4].map((x) => new Line(v(x, 0, 0), v(x, 2, 0)));
    const g = new GordonSurface(uCurves, vCurves);
    expect(g.parametersU[1]).toBeCloseTo(0.25, 9);
    expect(g.parametersV[1]).toBeCloseTo(0.5, 9);
    expectPoint(g.pointAt(0.25, 0.5), v(1, 1, 0), 1e-9);
  });

  it('refuses curves that cross, but at different normalized parameters', () => {
    const skewed = (j: number) => new ParametricCurve((t) => v(2 * t ** (1 + j), j, 0), 0, 1);
    const uCurves = [0, 1, 2].map(skewed);
    const vCurves = [0, 1, 2].map((x) => new Line(v(x, 0, 0), v(x, 2, 0)));
    expect(() => new GordonSurface(uCurves, vCurves)).toThrow(RangeError);
  });

  it('rejects a network that does not cross', () => {
    const a = [new Line(v(0, 0), v(1, 0)), new Line(v(0, 1), v(1, 1))];
    const b = [new Line(v(0, 0), v(0, 1)), new Line(v(5, 0), v(5, 1))];
    expect(() => new GordonSurface(a, b)).toThrow(RangeError);
    expect(() => new GordonSurface([a[0]!], b)).toThrow(RangeError);
    expect(() => new GordonSurface(a, [new Line(v(0, 0), v(0, 1)), new Line(v(1, 0), v(1, 1))], { parametersU: [0] })).toThrow(RangeError);
    expect(() => new GordonSurface(a, [new Line(v(0, 0), v(0, 1)), new Line(v(1, 0), v(1, 1))], { parametersU: [0.1, 1] })).toThrow(RangeError);
  });

  it('has a sweep along a straight path equal to an extrusion', () => {
    const profile = new BezierCurve([v(0, 0, 0), v(1, 0, 1), v(2, 0, 0)]);
    const path = new Line(v(0, 0, 0), v(0, 5, 0));
    const s = new SweepSurface(profile, path);
    const e = new ExtrusionSurface(profile, v(0, 5, 0));
    expect(s.equals(e, 1e-9)).toBe(true);
    expect(s.area / e.area).toBeCloseTo(1, 7);
    const t = new SweepSurface(profile, path, { mode: 'translation' });
    expect(t.equals(e, 1e-12)).toBe(true);
    expect(t.derivativeAt(0.5, 0.5, 1, 1).length).toBe(0);
  });

  it('turns a sweep with the path, but not a translation', () => {
    const profile = new Line(v(1, 0, 0), v(2, 0, 0));
    const path = new CircularArc(Frame3.WORLD, 1, 0, PI / 2);
    const turning = new SweepSurface(new Line(v(1, 0, 0), v(1, 0, 1)), path);
    // The profile stands on the circle, so it stays vertical all the way round.
    expectPoint(turning.pointAt(1, PI / 2), v(0, 1, 1), 1e-8);
    const moving = new SweepSurface(profile, path, { mode: 'translation' });
    expectPoint(moving.pointAt(1, PI / 2), v(1, 1, 0), 1e-12);
  });

  it('has a pipe with the area of a torus, on a closed path', () => {
    const [R, r] = [5, 1];
    const pipe = new PipeSurface(CircularArc.circle(Vector3.ZERO, R), r);
    expect(pipe.area).toBeCloseTo(4 * PI * PI * R * r, 6);
    expect(pipe.isPeriodicU).toBe(true);
    expect(pipe.isClosedV).toBe(true);
    expect(pipe.enclosesVolume).toBe(true);
    expect(pipe.enclosedVolume()).toBeCloseTo(2 * PI * PI * R * r * r, 5);
    const n = pipe.normalAt(0.3, 0.4);
    const axis = CircularArc.circle(Vector3.ZERO, R).pointAt(0.4);
    expectPoint(n, pipe.pointAt(0.3, 0.4).sub(axis).normalized(), 1e-6);
    expect(() => new PipeSurface(new Line(v(0, 0), v(1, 0)), 0)).toThrow(RangeError);
  });

  it('keeps every section of a pipe at the radius, and perpendicular to the path', () => {
    const path = BSplineCurve.interpolate([v(0, 0), v(5, 3), v(9, 0, 4), v(14, 3, 2)]);
    const pipe = new PipeSurface(path, 0.5);
    for (const b of [0, 0.2, 0.61, 1]) {
      for (const a of [0, 1.1, 4]) {
        const offset = pipe.pointAt(a, b).sub(path.pointAt(b));
        expect(offset.length).toBeCloseTo(0.5, 12);
        expect(offset.dot(path.tangentAt(b))).toBeCloseTo(0, 9);
      }
    }
  });

  it('has rotation-minimizing frames that do not twist around a planar curve', () => {
    const frames = new RotationMinimizingFrames(new EllipticalArc(Frame3.WORLD, 3, 2, 0, 3));
    for (const t of [0, 1, 2.5, 3]) {
      const f = frames.frameAt(t);
      expect(f.tangent.dot(f.normal)).toBeCloseTo(0, 9);
      expect(Math.abs(f.binormal.z)).toBeCloseTo(1, 8);
    }
    const straight = new RotationMinimizingFrames(new Polyline([v(0, 0), v(1, 0), v(1, 1)]));
    expect(straight.frameAt(1.5).tangent.y).toBeCloseTo(1, 12);
  });
});

describe('CompositeSurface', () => {
  const a = new PlaneSurface(Frame3.WORLD, { u: [0, 2], v: [0, 1] });
  const b = new PlaneSurface(Frame3.WORLD, { u: [2, 3], v: [0, 1] });
  const composite = new CompositeSurface([[a], [b]]);

  it('adds up the domains, and maps parameters to the pieces', () => {
    expect(composite.domainU.end).toBe(3);
    expect(composite.domainV.end).toBe(1);
    expectPoint(composite.pointAt(2.5, 0.5), v(2.5, 0.5));
    expect(composite.area).toBeCloseTo(3, 12);
    expect(composite.breakpointsU()).toEqual([0, 2, 3]);
    expect(composite.derivativeAt(2.5, 0.5, 1, 0).x).toBe(1);
  });

  it('keeps its class when transformed', () => {
    const moved = composite.transformed(Transform3.translate(v(1, 0, 0)));
    expect(moved).toBeInstanceOf(CompositeSurface);
    expectPoint(moved.pointAt(2.5, 0.5), v(3.5, 0.5));
  });

  it('rejects gaps and mismatches', () => {
    const far = new PlaneSurface(Frame3.WORLD, { u: [5, 6], v: [0, 1] });
    expect(() => new CompositeSurface([[a], [far]])).toThrow(RangeError);
    expect(() => new CompositeSurface([])).toThrow(RangeError);
    expect(() => new CompositeSurface([[a], []])).toThrow(RangeError);
    const shorter = new PlaneSurface(Frame3.WORLD, { u: [2, 3], v: [0, 0.5] });
    expect(() => new CompositeSurface([[a], [shorter]])).toThrow(RangeError);
    const top = new PlaneSurface(Frame3.WORLD, { u: [0, 1], v: [1, 2] });
    expect(() => new CompositeSurface([[a, top]])).toThrow(RangeError);
    const lifted = new PlaneSurface(Frame3.fromNormal(v(0, 0, 5)), { u: [0, 2], v: [1, 2] });
    expect(() => new CompositeSurface([[a, lifted]])).toThrow(RangeError);
  });
});

describe('ParametricSurface', () => {
  it('uses a derivative function if given, otherwise numeric ones', () => {
    const fn = (a: number, b: number) => v(a, b, a * a * b);
    const numeric = new ParametricSurface(fn, 0, 2, 0, 2);
    expectPoint(numeric.derivativeAt(1, 1, 1, 0), v(1, 0, 2), 1e-8);
    expectPoint(numeric.derivativeAt(1, 1, 2, 1), v(0, 0, 2), 1e-4);
    expect(() => numeric.derivativeAt(1, 1, 5, 0)).toThrow(RangeError);
    const exact = new ParametricSurface(fn, 0, 2, 0, 2, { derivative: () => v(9, 9, 9) });
    expectPoint(exact.derivativeAt(1, 1, 1, 0), v(9, 9, 9));
    expect(() => new ParametricSurface(fn, 1, 1, 0, 1)).toThrow(RangeError);
  });

  it('honours breakpoints and periodic directions', () => {
    const tent = (a: number, b: number) => v(a, b, Math.abs(a - 1));
    const s = new ParametricSurface(tent, 0, 2, 0, 1, { breakpointsU: [1, 5] });
    expect(s.breakpointsU()).toEqual([0, 1, 2]);
    expectPoint(s.derivativeAt(1, 0.5, 1, 0), v(1, 0, 1), 1e-7);
    expectPoint(s.derivativeAt(0.5, 0.5, 1, 0), v(1, 0, -1), 1e-7);
    expect(s.area).toBeCloseTo(2 * Math.SQRT2, 9);
    const nurbs = s.toNurbs(1e-6);
    expect(nurbs.distanceTo(v(1, 0.5, 0))).toBeLessThan(1e-6);
    const wrapped = new ParametricSurface((a, b) => v(Math.cos(a), Math.sin(a), b), 0, 2 * PI, 0, 1, { periodicU: true });
    expectPoint(wrapped.pointAt(7, 0.5), wrapped.pointAt(7 - 2 * PI, 0.5), 1e-12);
    expect(wrapped.isoCurveU(1).isPeriodic).toBe(false);
    expect(wrapped.isoCurveV(0.5).isPeriodic).toBe(true);
  });
});

describe('general forms', () => {
  const base = new TorusSurface(tilted, 3, 1);

  it('flattens nested trims, and clips breakpoints', () => {
    const loft = new LoftSurface([new Line(v(0, 0), v(1, 0)), new Line(v(0, 1), v(1, 1, 1)), new Line(v(0, 2), v(1, 2))], { parameters: [0, 1, 2] });
    const first = loft.trim(0, 1, 0.5, 2) as TrimmedSurface;
    expect(first).toBeInstanceOf(TrimmedSurface);
    expect(first.breakpointsV()).toEqual([0.5, 1, 2]);
    const second = first.trim(0, 1, 1.5, 1.7) as TrimmedSurface;
    expect(second.surface).toBe(loft);
    expect(() => first.trim(0, 1, 1.9, 1.9)).toThrow(RangeError);
  });

  it('undoes a double reversal, and mirrors breakpoints', () => {
    const r = new ReversedSurface(base);
    expect(r.reversed()).toBe(base);
    const loft = new LoftSurface([new Line(v(0, 0), v(1, 0)), new Line(v(0, 1), v(1, 1, 1))]);
    const polyU = new ExtrusionSurface(new Polyline([v(0, 0), v(1, 0), v(1, 1), v(3, 1)]), Vector3.Z);
    expect(polyU.reversed().breakpointsU()).toEqual([0, 1, 2, 3]);
    expect(loft.reversed().degreeU).toBeUndefined();
    const nurbs = BSplineSurface.uniform(2, 2, Array.from({ length: 4 }, (_, i) => Array.from({ length: 3 }, (_, j) => v(i, j, i * j))));
    const flipped = nurbs.reversed() as NurbsSurface;
    expect(flipped).toBeInstanceOf(NurbsSurface);
    expectPoint(flipped.pointAt(0.3, 0.4), nurbs.pointAt(0.7, 0.4), 1e-12);
  });

  it('composes transforms, and keeps NURBS exact under them', () => {
    const shear = new Transform3([1, 0.5, 0, 0, 1, 0, 0, 0, 1], Vector3.ZERO);
    const once = base.transformed(shear);
    expect(once).toBeInstanceOf(TransformedSurface);
    const twice = once.transformed(Transform3.translate(v(1, 1, 1))) as TransformedSurface;
    expect(twice.surface).toBe(base);
    const nurbs = once.toNurbs();
    expect(nurbs.isRational).toBe(true);
    expect(once.distanceTo(nurbs.pointAt(1.3, 2.1))).toBeLessThan(1e-7);
  });

  it('has an offset with the same parameters as its base', () => {
    const para = new EllipticParaboloidSurface(Frame3.WORLD, 2, 1, 0.5, { v: [0.3, 1] });
    const shell = new OffsetSurface(para, 0.25);
    expectPoint(shell.pointAt(1, 0.5), para.pointAt(1, 0.5).addScaled(para.normalAt(1, 0.5), 0.25));
    expect(shell.normalAt(1, 0.5).dot(para.normalAt(1, 0.5))).toBeCloseTo(1, 6);
  });

  it('has an iso-curve that behaves as a curve', () => {
    const iso = base.isoCurveU(0.7) as IsoCurve;
    expect(iso).toBeInstanceOf(IsoCurve);
    expect(iso.length).toBeCloseTo(2 * PI, 8);
    expect(iso.isPeriodic).toBe(true);
    expect(iso.curvatureAt(1)).toBeCloseTo(1, 7);
    const copy = curveFromJSON(JSON.parse(JSON.stringify(iso.toJSON())));
    expect(copy.equals(iso, 1e-9)).toBe(true);
    const plain = curveFromJSON(new Line(v(0, 0), v(1, 0)).toJSON());
    expect(plain).toBeInstanceOf(Line);
  });

  it('negates curvature when reversed, and keeps the Gaussian curvature', () => {
    const ball = new EllipsoidSurface(tilted, 3, 2, 1.5);
    const flipped = ball.reversed();
    const [a, b] = [1.2, 0.4];
    const mirrored = ball.domainU.start + ball.domainU.end - a;
    expect(flipped.meanCurvatureAt(mirrored, b)).toBeCloseTo(-ball.meanCurvatureAt(a, b), 9);
    expect(flipped.gaussianCurvatureAt(mirrored, b)).toBeCloseTo(ball.gaussianCurvatureAt(a, b), 9);
  });

  it('uses tangentU as the first direction where all directions are principal', () => {
    const sphere = new SphereSurface(tilted, 2);
    expectPoint(sphere.principalCurvaturesAt(1, 0.5).direction1, sphere.tangentUAt(1, 0.5), 1e-9);
  });

  it('refuses unknown JSON', () => {
    expect(() => surfaceFromJSON({ type: 'Nonsense' })).toThrow(TypeError);
  });
});

describe('CardinalSplines', () => {
  it('is 1 at its own knot and 0 at the others', () => {
    const c = new CardinalSplines([0, 1, 3, 4]);
    expect(c.count).toBe(4);
    [0, 1, 3, 4].forEach((k, i) => {
      c.weights(k).forEach((w, j) => expect(w).toBeCloseTo(i === j ? 1 : 0, 12));
    });
  });

  it('sums to one, so constants are reproduced', () => {
    const c = new CardinalSplines([0, 1, 3, 4]);
    for (const t of [0.5, 1.5, 2.2, 3.9]) {
      expect(c.weights(t).reduce((a, b) => a + b, 0)).toBeCloseTo(1, 12);
      expect(c.weights(t, 1).reduce((a, b) => a + b, 0)).toBeCloseTo(0, 12);
    }
  });

  it('is natural at the ends, and periodic if asked', () => {
    const c = new CardinalSplines([0, 1, 3, 4]);
    c.weights(0, 2).forEach((w) => expect(w).toBeCloseTo(0, 12));
    const p = new CardinalSplines([0, 1, 2, 4], true);
    expect(p.count).toBe(3);
    const a = p.weights(0, 1);
    const b = p.weights(4, 1);
    a.forEach((w, i) => expect(w).toBeCloseTo(b[i]!, 10));
  });

  it('rejects bad knots', () => {
    expect(() => new CardinalSplines([0])).toThrow(RangeError);
    expect(() => new CardinalSplines([0, 0, 1])).toThrow(RangeError);
    expect(() => new CardinalSplines([0, 1], true)).toThrow(RangeError);
  });
});
