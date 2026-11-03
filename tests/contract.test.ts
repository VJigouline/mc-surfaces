import { describe, expect, it } from 'vitest';
import { Transform3, Vector3, surfaceFromJSON, type ISurface, type SurfaceJSON } from '../src/index.js';
import { catalog, notSerializable } from './catalog.js';
import { bruteArea, expectPoint, smoothParameters, v } from './helpers.js';

/** Central differences with a span-aware step, independent of the code under test. */
function finite(s: ISurface, a: number, b: number, orderU: number, orderV: number): Vector3 {
  const hu = (s.domainU.end - s.domainU.start) * 1e-4;
  const hv = (s.domainV.end - s.domainV.start) * 1e-4;
  const p = (x: number, y: number) => s.pointAt(x, y);
  if (orderU === 1 && orderV === 0) {
    return p(a + hu, b).sub(p(a - hu, b)).scale(1 / (2 * hu));
  }
  if (orderU === 0 && orderV === 1) {
    return p(a, b + hv).sub(p(a, b - hv)).scale(1 / (2 * hv));
  }
  if (orderU === 2) {
    return p(a + hu, b).add(p(a - hu, b)).sub(p(a, b).scale(2)).scale(1 / (hu * hu));
  }
  if (orderV === 2) {
    return p(a, b + hv).add(p(a, b - hv)).sub(p(a, b).scale(2)).scale(1 / (hv * hv));
  }
  return p(a + hu, b + hv).sub(p(a + hu, b - hv)).sub(p(a - hu, b + hv)).add(p(a - hu, b - hv)).scale(1 / (4 * hu * hv));
}

const rotateAndMove = Transform3.translate(v(5, -3, 2)).compose(Transform3.rotate(v(1, 2, -1), 0.8));
const span = (d: { start: number; end: number }, s: number) => d.start + (d.end - d.start) * s;

describe.each(Object.entries(catalog))('%s', (name, make) => {
  const surface = make();
  const size = Math.max(surface.boundingBox().diagonal, 1);
  const { domainU: du, domainV: dv } = surface;

  it('has valid domains and breakpoints', () => {
    expect(du.end).toBeGreaterThan(du.start);
    expect(dv.end).toBeGreaterThan(dv.start);
    const [bu, bv] = [surface.breakpointsU(), surface.breakpointsV()];
    expect(bu[0]).toBe(du.start);
    expect(bu[bu.length - 1]).toBe(du.end);
    expect(bv[0]).toBe(dv.start);
    expect(bv[bv.length - 1]).toBe(dv.end);
  });

  it('clamps or wraps parameters outside the domain', () => {
    const b = span(dv, 0.4);
    const below = surface.pointAt(du.start - (du.end - du.start) / 4, b);
    const expected = surface.isPeriodicU ? surface.pointAt(span(du, 0.75), b) : surface.pointAt(du.start, b);
    expectPoint(below, expected, 1e-8 * size);
    expect(() => surface.pointAt(NaN, b)).toThrow(RangeError);
  });

  it('has derivatives matching finite differences', () => {
    for (const [a, b] of smoothParameters(surface, 4)) {
      for (const [ou, ov, tol] of [[1, 0, 1e-5], [0, 1, 1e-5], [2, 0, 1e-3], [0, 2, 1e-3], [1, 1, 1e-3]] as const) {
        const d = surface.derivativeAt(a, b, ou, ov);
        expectPoint(d, finite(surface, a, b, ou, ov), tol * Math.max(d.length, 1));
      }
    }
  });

  it('rejects invalid derivative orders', () => {
    const [a, b] = [span(du, 0.5), span(dv, 0.5)];
    expect(() => surface.derivativeAt(a, b, -1, 0)).toThrow(RangeError);
    expect(() => surface.derivativeAt(a, b, 1.5, 0)).toThrow(RangeError);
    expectPoint(surface.derivativeAt(a, b, 0, 0), surface.pointAt(a, b));
  });

  it('has a unit normal along Su x Sv', () => {
    for (const [a, b] of smoothParameters(surface, 4)) {
      const n = surface.normalAt(a, b);
      expect(n.length).toBeCloseTo(1, 9);
      const c = surface.derivativeAt(a, b, 1, 0).cross(surface.derivativeAt(a, b, 0, 1));
      if (c.length > 1e-6) {
        expectPoint(n, c.normalized(), 1e-6);
      }
      const f = surface.frameAt(a, b);
      expect(f.tangentU.dot(n)).toBeCloseTo(0, 6);
      expect(f.tangentV.dot(n)).toBeCloseTo(0, 6);
    }
  });

  it('has curvatures consistent with its fundamental forms', () => {
    for (const [a, b] of smoothParameters(surface, 3)) {
      const su = surface.derivativeAt(a, b, 1, 0);
      const sv = surface.derivativeAt(a, b, 0, 1);
      if (su.cross(sv).length < 1e-6) {
        continue;
      }
      const n = surface.normalAt(a, b);
      const { E, F, G } = surface.firstFundamentalFormAt(a, b);
      const { L, M, N } = surface.secondFundamentalFormAt(a, b);
      expect(E).toBeCloseTo(su.dot(su), 9);
      expect(F).toBeCloseTo(su.dot(sv), 9);
      expect(G).toBeCloseTo(sv.dot(sv), 9);
      const p = surface.principalCurvaturesAt(a, b);
      expect(p.k1).toBeGreaterThanOrEqual(p.k2);
      const det = E * G - F * F;
      const scale = Math.max(Math.abs(p.k1), Math.abs(p.k2), 1);
      expect(Math.abs(p.k1 * p.k2 - (L * N - M * M) / det)).toBeLessThan(1e-6 * scale * scale);
      expect(surface.gaussianCurvatureAt(a, b)).toBeCloseTo(p.k1 * p.k2, 6);
      expect(surface.meanCurvatureAt(a, b)).toBeCloseTo((p.k1 + p.k2) / 2, 6);
      expect(p.direction1.length).toBeCloseTo(1, 9);
      expect(p.direction1.dot(n)).toBeCloseTo(0, 6);
      expect(p.direction2.dot(p.direction1)).toBeCloseTo(0, 6);
      // Normal curvature along the first principal direction, written in (u, v) terms.
      if (p.k1 - p.k2 > 1e-6) {
        const [x, y] = [p.direction1.dot(su), p.direction1.dot(sv)];
        const c = (G * x - F * y) / det;
        const d = (E * y - F * x) / det;
        const normalCurvature = (L * c * c + 2 * M * c * d + N * d * d) / (E * c * c + 2 * F * c * d + G * d * d);
        expect(Math.abs(normalCurvature - p.k1)).toBeLessThan(1e-5 * scale);
      }
    }
  });

  it('has iso-curves and boundaries that lie on it', () => {
    const [a, b] = [span(du, 0.37), span(dv, 0.61)];
    const alongV = surface.isoCurveU(a);
    const alongU = surface.isoCurveV(b);
    for (const s of [0, 0.3, 0.8, 1]) {
      expectPoint(alongV.pointAt(span(alongV.domain, s)), surface.pointAt(a, span(dv, s)), 1e-8 * size);
      expectPoint(alongU.pointAt(span(alongU.domain, s)), surface.pointAt(span(du, s), b), 1e-8 * size);
    }
    const [e0, e1, e2, e3] = surface.boundaryCurves();
    expectPoint(e0.startPoint, surface.pointAt(du.start, dv.start), 1e-8 * size);
    expectPoint(e1.endPoint, surface.pointAt(du.end, dv.end), 1e-8 * size);
    expectPoint(e2.startPoint, surface.pointAt(du.start, dv.end), 1e-8 * size);
    expectPoint(e3.endPoint, surface.pointAt(du.start, dv.end), 1e-8 * size);
  });

  it('measures area like a fine mesh', { timeout: 60000 }, () => {
    const brute = bruteArea(surface, 300);
    expect(surface.area / brute).toBeCloseTo(1, 3);
    const half = surface.areaBetween(du.start, span(du, 0.5), dv.start, dv.end);
    const rest = surface.areaBetween(span(du, 0.5), du.end, dv.start, dv.end);
    expect((half + rest) / surface.area).toBeCloseTo(1, 7);
    expect(surface.areaBetween(span(du, 0.5), span(du, 0.5), dv.start, dv.end)).toBe(0);
  });

  it('has a bounding box that contains it tightly', () => {
    const box = surface.boundingBox();
    const dense = surface.sample(120, 120).flat();
    for (const p of dense) {
      expect(box.contains(p, 2e-5 * size)).toBe(true);
    }
    for (const axis of ['x', 'y', 'z'] as const) {
      const values = dense.map((p) => p[axis]);
      expect(Math.min(...values) - box.min[axis]).toBeLessThan(2e-3 * size);
      expect(box.max[axis] - Math.max(...values)).toBeLessThan(2e-3 * size);
    }
  });

  it('finds closest points', { timeout: 60000 }, () => {
    for (const [a, b] of [[0.31, 0.27], [0.77, 0.58], [0.5, 0.91]] as const) {
      const on = surface.pointAt(span(du, a), span(dv, b));
      expect(surface.distanceTo(on)).toBeLessThan(1e-5 * size);
      const out = on.add(v(0.3, -0.7, 0.5).scale(size * 0.4));
      const found = surface.closestPoint(out);
      expectPoint(found.point, surface.pointAt(found.u, found.v), 1e-9 * size);
      expect(found.distance).toBeCloseTo(found.point.distanceTo(out), 9);
      // No sample may beat the reported distance.
      for (const p of surface.sample(60, 60).flat()) {
        expect(found.distance).toBeLessThanOrEqual(p.distanceTo(out) + 1e-6 * size);
      }
    }
  });

  it('tessellates within tolerance', { timeout: 60000 }, () => {
    const chordTolerance = size / 100;
    const mesh = surface.tessellate({ chordTolerance });
    expect(mesh.triangles.length % 3).toBe(0);
    expect(mesh.triangles.length).toBeGreaterThan(0);
    expect(mesh.normals.length).toBe(mesh.points.length);
    const stride = Math.max(1, Math.floor(mesh.triangles.length / 3 / 500));
    for (let t = 0; t < mesh.triangles.length / 3; t += stride) {
      const [ia, ib, ic] = [mesh.triangles[3 * t]!, mesh.triangles[3 * t + 1]!, mesh.triangles[3 * t + 2]!];
      const [a, b, c] = [mesh.points[ia]!, mesh.points[ib]!, mesh.points[ic]!];
      const [pa, pb, pc] = [mesh.parameters[ia]!, mesh.parameters[ib]!, mesh.parameters[ic]!];
      expectPoint(a, surface.pointAt(pa.u, pa.v), 1e-9 * size);
      const centre = a.add(b).add(c).scale(1 / 3);
      const mid = surface.pointAt((pa.u + pb.u + pc.u) / 3, (pa.v + pb.v + pc.v) / 3);
      expect(centre.distanceTo(mid)).toBeLessThan(1.5 * chordTolerance);
      // Winding agrees with the normal.
      const facing = b.sub(a).cross(c.sub(a));
      if (facing.length > 1e-9 * size * size) {
        const pointing = mesh.normals[ia]!.add(mesh.normals[ib]!).add(mesh.normals[ic]!);
        expect(facing.dot(pointing)).toBeGreaterThan(0);
      }
    }
    expect(surface.tessellate({ chordTolerance: size / 1000 }).points.length).toBeGreaterThanOrEqual(mesh.points.length);
    expect(() => surface.tessellate({ chordTolerance: 0 })).toThrow(RangeError);
  });

  it('samples a grid', () => {
    const grid = surface.sample(3, 4);
    expect(grid.length).toBe(3);
    expect(grid[0]!.length).toBe(4);
    expectPoint(grid[2]![3]!, surface.pointAt(du.end, dv.end), 1e-9 * size);
    expect(() => surface.sample(1, 4)).toThrow(RangeError);
  });

  it('keeps parameters and geometry when trimmed', { timeout: 60000 }, () => {
    const [a, b, c, d] = [span(du, 0.2), span(du, 0.7), span(dv, 0.3), span(dv, 0.9)];
    const piece = surface.trim(a, b, c, d);
    expect(piece.domainU.start).toBeCloseTo(a, 12);
    expect(piece.domainV.end).toBeCloseTo(d, 12);
    for (const [x, y] of [[0.3, 0.4], [0.5, 0.8], [0.65, 0.35]] as const) {
      expectPoint(piece.pointAt(span(du, x), span(dv, y)), surface.pointAt(span(du, x), span(dv, y)), 1e-7 * size);
    }
    expect(() => surface.trim(a, a, c, d)).toThrow(RangeError);
    const [before, after] = surface.splitU(span(du, 0.4));
    expectPoint(before.pointAt(span(du, 0.2), span(dv, 0.5)), surface.pointAt(span(du, 0.2), span(dv, 0.5)), 1e-7 * size);
    expectPoint(after.pointAt(span(du, 0.8), span(dv, 0.5)), surface.pointAt(span(du, 0.8), span(dv, 0.5)), 1e-7 * size);
    const [low, high] = surface.splitV(span(dv, 0.6));
    expect((low.area + high.area) / surface.area).toBeCloseTo(1, 5);
    expect(() => surface.splitU(du.start)).toThrow(RangeError);
    expect(() => surface.splitV(dv.end)).toThrow(RangeError);
  });

  it('flips its normal when reversed', { timeout: 60000 }, () => {
    const flipped = surface.reversed();
    const [a, b] = [span(du, 0.37), span(dv, 0.43)];
    const mirrored = du.start + du.end - a;
    expectPoint(flipped.pointAt(mirrored, b), surface.pointAt(a, b), 1e-9 * size);
    expectPoint(flipped.normalAt(mirrored, b), surface.normalAt(a, b).negate(), 1e-7);
    expect(flipped.area / surface.area).toBeCloseTo(1, 5);
  });

  it('keeps parameters under rigid, scaling and shearing transforms', () => {
    const shear = new Transform3([1, 0.3, 0, 0, 1.2, 0, 0.1, 0, 0.8], v(1, 0, 0));
    for (const t of [rotateAndMove, Transform3.scale(2.5), Transform3.mirror(v(1, 1, 0)), shear]) {
      const moved = surface.transformed(t);
      const again = moved.transformed(rotateAndMove);
      for (const [x, y] of [[0.31, 0.27], [0.77, 0.58]] as const) {
        const [a, b] = [span(du, x), span(dv, y)];
        expectPoint(moved.pointAt(a, b), t.applyToPoint(surface.pointAt(a, b)), 1e-6 * size);
        expectPoint(again.pointAt(a, b), rotateAndMove.applyToPoint(t.applyToPoint(surface.pointAt(a, b))), 1e-6 * size);
        const d = moved.derivativeAt(a, b, 1, 0);
        expectPoint(d, t.applyToVector(surface.derivativeAt(a, b, 1, 0)), 1e-5 * Math.max(d.length, 1));
      }
    }
  });

  it('offsets along the normal', { timeout: 60000 }, () => {
    const distance = 0.02 * size;
    const parallel = surface.offset(distance);
    for (const [x, y] of [[0.31, 0.27], [0.77, 0.58]] as const) {
      const [a, b] = [span(du, x), span(dv, y)];
      const target = surface.pointAt(a, b).addScaled(surface.normalAt(a, b), distance);
      expect(parallel.distanceTo(target)).toBeLessThan(1e-4 * size);
    }
  });

  it('converts to NURBS within tolerance', { timeout: 120000 }, () => {
    const tolerance = 1e-3 * size;
    const nurbs = surface.toNurbs(tolerance);
    expect(nurbs.controlPoints.length).toBeGreaterThan(1);
    for (const [x, y] of [[0.13, 0.27], [0.5, 0.5], [0.77, 0.58], [0.9, 0.1]] as const) {
      expect(nurbs.distanceTo(surface.pointAt(span(du, x), span(dv, y)))).toBeLessThan(3 * tolerance);
    }
    for (const [x, y] of [[0.3, 0.7], [0.6, 0.2]] as const) {
      const on = nurbs.pointAt(span(nurbs.domainU, x), span(nurbs.domainV, y));
      expect(surface.distanceTo(on)).toBeLessThan(3 * tolerance);
    }
  });

  it('has a plane section that lies on it and on the plane', { timeout: 60000 }, () => {
    const mid = surface.pointAt(span(du, 0.53), span(dv, 0.47));
    const normal = v(0.3, 0.5, 0.8).normalized();
    const sections = surface.intersectPlane(mid, normal);
    for (const section of sections) {
      expect(section.points.length).toBe(section.parameters.length);
      section.points.forEach((p, i) => {
        expect(Math.abs(p.sub(mid).dot(normal))).toBeLessThan(1e-7 * size);
        const q = section.parameters[i]!;
        expectPoint(p, surface.pointAt(q.u, q.v), 1e-9 * size);
      });
    }
    expect(sections.some((s) => s.points.some((p) => p.distanceTo(mid) < size / 10))).toBe(true);
  });

  it('intersects a line through a known point', { timeout: 60000 }, () => {
    const [a, b] = [span(du, 0.53), span(dv, 0.47)];
    const target = surface.pointAt(a, b);
    const n = surface.normalAt(a, b);
    const hits = surface.intersectLine(target.addScaled(n, size), n, 1e-6 * size);
    expect(hits.some((h) => h.point.distanceTo(target) < 1e-4 * size)).toBe(true);
    for (const h of hits) {
      expectPoint(h.point, surface.pointAt(h.u, h.v), 1e-5 * size);
    }
  });

  it('reports consistent closure', () => {
    const b = span(dv, 0.5);
    if (surface.isClosedU) {
      expectPoint(surface.pointAt(du.start, b), surface.pointAt(du.end, b), 1e-7 * size);
    }
    if (surface.isPeriodicU) {
      expect(surface.isClosedU).toBe(true);
    }
    if (surface.isPeriodicV) {
      expect(surface.isClosedV).toBe(true);
    }
    expect(surface.enclosedVolume() === undefined).toBe(!surface.enclosesVolume);
  });

  it('equals itself and not a moved copy', () => {
    expect(surface.equals(surface)).toBe(true);
    expect(surface.equals(surface.transformed(Transform3.translate(v(size, 0, 0))))).toBe(false);
  });

  it(notSerializable.has(name) ? 'refuses JSON' : 'round-trips through JSON', () => {
    if (notSerializable.has(name)) {
      expect(() => surface.toJSON()).toThrow(TypeError);
      return;
    }
    const json = JSON.parse(JSON.stringify(surface.toJSON())) as SurfaceJSON;
    const copy = surfaceFromJSON(json);
    expect(copy.type).toBe(surface.type);
    expect(surface.equals(copy, 1e-7 * size)).toBe(true);
    expect(copy.domainU.end).toBeCloseTo(du.end, 12);
  });

  it('has a consistent centroid', { timeout: 60000 }, () => {
    const c = surface.centroid();
    expect(surface.boundingBox().contains(c, 1e-6 * size)).toBe(true);
    expectPoint(surface.transformed(Transform3.translate(v(3, 4, 5))).centroid(), c.add(v(3, 4, 5)), 1e-5 * size);
  });
});
