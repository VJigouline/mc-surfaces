import { describe, expect, it } from 'vitest';
import { Transform3, Vector3, SphereSurface, type ISurface } from '../src/index.js';
import { catalog } from './catalog.js';
import { expectPoint } from './helpers.js';

/** Model sizes the documentation says are supported. */
const SCALES = [1e-3, 1e6];

describe.each(SCALES)('at model scale %d', (scale) => {
  const scaleUp = Transform3.scale(scale);

  describe.each(Object.entries(catalog))('%s', (_name, make) => {
    const base = make();
    const surface: ISurface = base.transformed(scaleUp);
    const { domainU: du, domainV: dv } = surface;
    const a = du.start + (du.end - du.start) * 0.37;
    const b = dv.start + (dv.end - dv.start) * 0.41;

    it('scales its measures, promptly', { timeout: 15000 }, () => {
      expect(surface.area / (base.area * scale * scale)).toBeCloseTo(1, 5);
      expect(surface.boundingBox().diagonal / (base.boundingBox().diagonal * scale)).toBeCloseTo(1, 3);
      expect(surface.isClosedU).toBe(base.isClosedU);
      expect(surface.isClosedV).toBe(base.isClosedV);
      expect(surface.enclosesVolume).toBe(base.enclosesVolume);
    });

    it('keeps its normal and scales its curvature', { timeout: 15000 }, () => {
      expectPoint(surface.normalAt(a, b), base.normalAt(a, b), 1e-5);
      const small = base.principalCurvaturesAt(a, b);
      const big = surface.principalCurvaturesAt(a, b);
      expect(Math.abs(big.k1 * scale - small.k1)).toBeLessThan(1e-3 * Math.max(1, Math.abs(small.k1)));
    });

    it('finds closest points, and tessellates in proportion', { timeout: 15000 }, () => {
      expect(surface.distanceTo(surface.pointAt(a, b)) / scale).toBeLessThan(1e-5);
      const size = base.boundingBox().diagonal;
      const coarse = base.tessellate({ chordTolerance: size / 100 }).points.length;
      const scaled = surface.tessellate({ chordTolerance: (size * scale) / 100 }).points.length;
      expect(scaled).toBe(coarse);
    });
  });

  it('intersects, given tolerances in proportion', { timeout: 15000 }, () => {
    const sphere = SphereSurface.centredAt(Vector3.ZERO, scale);
    const hits = sphere.intersectLine(new Vector3(-2 * scale, 0, 0), Vector3.X, 1e-7 * scale);
    expect(hits).toHaveLength(2);
    expectPoint(hits[1]!.point, new Vector3(scale, 0, 0), 1e-5 * scale);
    const other = SphereSurface.centredAt(new Vector3(0, 0, scale), scale);
    const curves = sphere.intersect(other, { tolerance: 1e-7 * scale });
    expect(curves).toHaveLength(1);
    expect(curves[0]!.closed).toBe(true);
    expect(sphere.intersectPlane(Vector3.ZERO, Vector3.Z)).toHaveLength(1);
  });
});
