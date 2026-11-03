import { expect } from 'vitest';
import { Vector3, type ISurface } from '../src/index.js';

export const v = (x: number, y: number, z = 0) => new Vector3(x, y, z);

export function expectPoint(actual: Vector3, expected: Vector3, tolerance = 1e-9): void {
  const distance = actual.distanceTo(expected);
  if (!(distance <= tolerance)) {
    expect.fail(`(${actual.toArray()}) is ${distance} from (${expected.toArray()}), tolerance ${tolerance}`);
  }
}

/** Evenly spaced parameters over a domain, strictly inside it. */
export function inner(domain: { start: number; end: number }, count: number): number[] {
  return Array.from({ length: count }, (_, i) => domain.start + ((domain.end - domain.start) * (i + 0.5)) / count);
}

/** Parameter pairs a small margin away from every breakpoint, for tests of smooth behaviour. */
export function smoothParameters(surface: ISurface, count: number): [number, number][] {
  const breaksU = surface.breakpointsU();
  const breaksV = surface.breakpointsV();
  const marginU = (surface.domainU.end - surface.domainU.start) * 4e-3;
  const marginV = (surface.domainV.end - surface.domainV.start) * 4e-3;
  const us = inner(surface.domainU, count).filter((u) => breaksU.every((b) => Math.abs(b - u) > marginU));
  const vs = inner(surface.domainV, count).filter((w) => breaksV.every((b) => Math.abs(b - w) > marginV));
  return us.flatMap((u) => vs.map((w) => [u, w] as [number, number]));
}

/** Area of a fine mesh of the surface, as an independent check of the integrator. */
export function bruteArea(surface: ISurface, segments = 400): number {
  const grid = surface.sample(segments + 1, segments + 1);
  let total = 0;
  for (let i = 0; i < segments; i++) {
    for (let j = 0; j < segments; j++) {
      const [a, b, c, d] = [grid[i]![j]!, grid[i + 1]![j]!, grid[i + 1]![j + 1]!, grid[i]![j + 1]!];
      total += b.sub(a).cross(c.sub(a)).length / 2 + c.sub(a).cross(d.sub(a)).length / 2;
    }
  }
  return total;
}
