import { CircularArc, Helix, Line } from 'mc-curves';
import { describe, expect, it } from 'vitest';
import {
  BSplineSurface,
  CylinderSurface,
  Frame3,
  PlaneSurface,
  SphereSurface,
  TorusSurface,
  Vector3,
} from '../src/index.js';
import { expectPoint, v } from './helpers.js';

const PI = Math.PI;
const sphere = SphereSurface.centredAt(Vector3.ZERO, 2);
const cylinder = new CylinderSurface(Frame3.fromNormal(v(0, 0, -5)), 1, { v: [0, 10] });

describe('intersectLine', () => {
  it('finds both crossings of a sphere, sorted along the line, in distance from its origin', () => {
    const hits = sphere.intersectLine(v(-5, 0, 0), Vector3.X);
    expect(hits).toHaveLength(2);
    expectPoint(hits[0]!.point, v(-2, 0, 0), 1e-7);
    expectPoint(hits[1]!.point, v(2, 0, 0), 1e-7);
    expect(hits[0]!.t).toBeCloseTo(3, 7);
    expect(hits[1]!.t).toBeCloseTo(7, 7);
    expectPoint(sphere.pointAt(hits[0]!.u, hits[0]!.v), hits[0]!.point, 1e-6);
  });

  it('misses a sphere it passes by', () => {
    expect(sphere.intersectLine(v(0, 3, 0), Vector3.X)).toHaveLength(0);
  });

  it('crosses through a seam and a pole', () => {
    expect(sphere.intersectLine(v(0, 0, -5), Vector3.Z)).toHaveLength(2);
    const hits = sphere.intersectLine(v(5, 0, 0), Vector3.X.negate());
    expect(hits).toHaveLength(2);
    expectPoint(hits[0]!.point, v(2, 0, 0), 1e-7);
  });

  it('finds a tangent touch in this clear case', () => {
    const hits = sphere.intersectLine(v(-5, 2, 0), Vector3.X);
    expect(hits.length).toBeGreaterThan(0);
    expect(hits[0]!.point.distanceTo(v(0, 2, 0))).toBeLessThan(1e-2);
  });

  it('crosses a plane only inside its rectangle', () => {
    const plate = new PlaneSurface(Frame3.WORLD, { u: [0, 1], v: [0, 1] });
    expect(plate.intersectLine(v(0.5, 0.5, 3), Vector3.Z)).toHaveLength(1);
    expect(plate.intersectLine(v(1.5, 0.5, 3), Vector3.Z)).toHaveLength(0);
  });

  it('rejects a zero direction', () => {
    expect(() => sphere.intersectLine(Vector3.ZERO, Vector3.ZERO)).toThrow(RangeError);
  });
});

describe('intersectCurve', () => {
  it('stops at the end of a finite curve', () => {
    const short = new Line(v(-5, 0, 0), v(0, 0, 0));
    const hits = sphere.intersectCurve(short);
    expect(hits).toHaveLength(1);
    expectPoint(hits[0]!.point, v(-2, 0, 0), 1e-7);
    expect(hits[0]!.t).toBeCloseTo(0.6, 7);
    expect(hits[0]!.distance).toBeLessThanOrEqual(1e-6);
  });

  it('skips a curve whose box is clear of the surface', () => {
    expect(sphere.intersectCurve(new Line(v(10, 10, 10), v(11, 10, 10)))).toHaveLength(0);
  });

  it('counts each turn of a helix through a plane', () => {
    const helix = new Helix(Frame3.WORLD, 2, 1, 5);
    const floor = PlaneSurface.centred(v(0, 0, 2.3), Vector3.Z, 10, 10);
    expect(floor.intersectCurve(helix)).toHaveLength(1);
    const wall = PlaneSurface.centred(v(0, 0, 2.5), Vector3.X, 10, 10, Vector3.Y);
    expect(wall.intersectCurve(helix)).toHaveLength(10);
  });

  it('crosses a circle with a cylinder at four points', () => {
    const ring = new CircularArc(Frame3.fromNormal(Vector3.ZERO, Vector3.Y, Vector3.X), 2, 0, 2 * PI);
    const hits = new CylinderSurface(Frame3.fromNormal(v(-5, 0, 0), Vector3.X), 1, { v: [0, 10] }).intersectCurve(ring);
    expect(hits).toHaveLength(4);
    for (const h of hits) {
      expect(Math.abs(h.point.z)).toBeCloseTo(1, 6);
      expect(Math.abs(h.point.x)).toBeCloseTo(Math.sqrt(3), 6);
    }
  });
});

describe('intersectPlane', () => {
  const circleRadius = (points: readonly Vector3[], centre: Vector3) => points.map((p) => p.distanceTo(centre));

  it('cuts a sphere in one closed circle', () => {
    const sections = sphere.intersectPlane(v(0, 0, 1), Vector3.Z);
    expect(sections).toHaveLength(1);
    expect(sections[0]!.closed).toBe(true);
    for (const r of circleRadius(sections[0]!.points, v(0, 0, 1))) {
      expect(r).toBeCloseTo(Math.sqrt(3), 8);
    }
  });

  it('cuts a sphere through its seam and poles', () => {
    const sections = sphere.intersectPlane(Vector3.ZERO, Vector3.Y);
    expect(sections.length).toBeGreaterThan(0);
    for (const s of sections) {
      for (const p of s.points) {
        expect(p.length).toBeCloseTo(2, 8);
        expect(Math.abs(p.y)).toBeLessThan(1e-8);
      }
    }
    const span = sections.flatMap((s) => s.points);
    expect(Math.max(...span.map((p) => p.x))).toBeCloseTo(2, 1);
    expect(Math.min(...span.map((p) => p.x))).toBeCloseTo(-2, 1);
  });

  it('misses a sphere it does not reach', () => {
    expect(sphere.intersectPlane(v(0, 0, 3), Vector3.Z)).toHaveLength(0);
  });

  it('cuts a torus through its axis in two circles, and through its middle in two more', () => {
    const torus = new TorusSurface(Frame3.WORLD, 3, 1);
    const across = torus.intersectPlane(Vector3.ZERO, Vector3.Y);
    expect(across).toHaveLength(2);
    for (const s of across) {
      expect(s.closed).toBe(true);
      const side = s.points[0]!.x > 0 ? 1 : -1;
      for (const r of circleRadius(s.points, v(3 * side, 0, 0))) {
        expect(r).toBeCloseTo(1, 8);
      }
    }
    const level = torus.intersectPlane(v(0, 0, 0.1), Vector3.Z);
    expect(level).toHaveLength(2);
  });

  it('cuts an open surface in an open polyline, with both ends on its edge', () => {
    const plate = new PlaneSurface(Frame3.WORLD, { u: [0, 4], v: [0, 2] });
    const sections = plate.intersectPlane(v(1, 0, 0), v(1, 0, 0.2));
    expect(sections).toHaveLength(1);
    expect(sections[0]!.closed).toBe(false);
    const ys = sections[0]!.points.map((p) => p.y);
    expect(Math.min(...ys)).toBeCloseTo(0, 9);
    expect(Math.max(...ys)).toBeCloseTo(2, 9);
  });

  it('cuts a cylinder at an angle in an ellipse', () => {
    const sections = cylinder.intersectPlane(Vector3.ZERO, v(0, 1, 1));
    expect(sections).toHaveLength(1);
    expect(sections[0]!.closed).toBe(true);
    for (const p of sections[0]!.points) {
      expect(Math.hypot(p.x, p.y)).toBeCloseTo(1, 8);
      expect(p.y + p.z).toBeCloseTo(0, 8);
    }
  });

  it('gives nothing for a surface lying in the plane, and rejects a zero normal', () => {
    const plate = new PlaneSurface(Frame3.WORLD, { u: [0, 4], v: [0, 2] });
    expect(plate.intersectPlane(Vector3.ZERO, Vector3.Z)).toHaveLength(0);
    expect(() => plate.intersectPlane(Vector3.ZERO, Vector3.ZERO)).toThrow(RangeError);
  });
});

describe('intersect', () => {
  it('cuts two spheres in a circle', () => {
    const other = SphereSurface.centredAt(v(0, 0, 2), 2);
    const curves = sphere.intersect(other);
    expect(curves).toHaveLength(1);
    expect(curves[0]!.closed).toBe(true);
    const centre = v(0, 0, 1);
    for (const p of curves[0]!.points) {
      expect(p.z).toBeCloseTo(1, 5);
      expect(p.distanceTo(centre)).toBeCloseTo(Math.sqrt(3), 5);
    }
    curves[0]!.parameters.forEach((q, i) => {
      expect(sphere.pointAt(q.u, q.v).distanceTo(curves[0]!.points[i]!)).toBeLessThan(2e-6);
      const r = curves[0]!.parametersOther[i]!;
      expect(other.pointAt(r.u, r.v).distanceTo(curves[0]!.points[i]!)).toBeLessThan(2e-6);
    });
  });

  it('cuts a plane through a sphere in a circle', () => {
    const plate = PlaneSurface.centred(v(0, 0, 0.5), Vector3.Z, 10, 10);
    const curves = plate.intersect(sphere);
    expect(curves).toHaveLength(1);
    for (const p of curves[0]!.points) {
      expect(p.z).toBeCloseTo(0.5, 6);
      expect(Math.hypot(p.x, p.y)).toBeCloseTo(Math.sqrt(3.75), 5);
    }
  });

  it('cuts two crossing cylinders in two closed loops', () => {
    const across = new CylinderSurface(Frame3.fromNormal(v(-5, 0, 0), Vector3.X), 0.5, { v: [0, 10] });
    const curves = cylinder.intersect(across);
    expect(curves).toHaveLength(2);
    for (const c of curves) {
      expect(c.closed).toBe(true);
      for (const p of c.points) {
        expect(Math.hypot(p.x, p.y)).toBeCloseTo(1, 5);
        expect(Math.hypot(p.y, p.z)).toBeCloseTo(0.5, 5);
      }
    }
  });

  it('returns an open curve across two patches', () => {
    const a = new PlaneSurface(Frame3.WORLD, { u: [0, 2], v: [0, 2] });
    const wall = PlaneSurface.centred(v(1, 1, 0), Vector3.X, 8, 8, Vector3.Y);
    const curves = a.intersect(wall);
    expect(curves).toHaveLength(1);
    expect(curves[0]!.closed).toBe(false);
    const ys = curves[0]!.points.map((p) => p.y);
    expect(Math.min(...ys)).toBeCloseTo(0, 6);
    expect(Math.max(...ys)).toBeCloseTo(2, 6);
  });

  it('cuts a free-form surface with a plane, onto both surfaces', () => {
    const net = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => v(i, j, Math.sin(i) * Math.cos(j))));
    const spline = BSplineSurface.uniform(3, 3, net);
    const plate = PlaneSurface.centred(v(2.5, 2.5, 0.1), Vector3.Z, 20, 20);
    const curves = spline.intersect(plate);
    expect(curves.length).toBeGreaterThan(0);
    for (const c of curves) {
      for (const q of c.parameters) {
        expect(spline.pointAt(q.u, q.v).z).toBeCloseTo(0.1, 5);
      }
    }
  });

  it('gives nothing for surfaces that do not meet, or do not cross', () => {
    expect(sphere.intersect(SphereSurface.centredAt(v(10, 0, 0), 1))).toHaveLength(0);
    expect(sphere.intersect(SphereSurface.centredAt(v(0, 0, 0), 1))).toHaveLength(0);
    const plate = new PlaneSurface(Frame3.WORLD, { u: [0, 1], v: [0, 1] });
    expect(plate.intersect(plate)).toHaveLength(0);
  });
});
