import type { Vector3 } from 'meercad-curves';
import { gridParameters } from '../math/numeric.js';
import type { ISurface, SurfaceMesh, SurfaceParameter, SurfaceTessellationOptions } from '../surfaces/ISurface.js';

const INITIAL_PIECES = 8;
const MAX_PASSES = 16;

function chordDeviation(p: Vector3, a: Vector3, b: Vector3): number {
  return p.distanceTo(a.lerp(b, 0.5));
}

function refined(values: readonly number[], split: ReadonlySet<number>): number[] {
  const result: number[] = [];
  values.forEach((value, i) => {
    result.push(value);
    if (split.has(i)) {
      result.push((value + values[i + 1]!) / 2);
    }
  });
  return result;
}

/**
 * The parameter lines of a tensor grid, refined until each cell's edges and diagonal are within
 * the tolerance. A whole line is added at a time, so the grid has no cracks to patch.
 */
function refineGrid(surface: ISurface, options: SurfaceTessellationOptions): { us: number[]; vs: number[] } {
  const tolerance = options.chordTolerance ?? 0.01;
  if (!(tolerance > 0)) {
    throw new RangeError('Chord tolerance must be positive.');
  }
  const maxEdge = options.maxEdgeLength ?? Infinity;
  const maxPoints = options.maxPoints ?? 1_000_000;
  let us = gridParameters(surface.breakpointsU(), INITIAL_PIECES);
  let vs = gridParameters(surface.breakpointsV(), INITIAL_PIECES);
  for (let pass = 0; pass < MAX_PASSES; pass++) {
    const grid = us.map((u) => vs.map((v) => surface.pointAt(u, v)));
    const splitU = new Set<number>();
    const splitV = new Set<number>();
    for (let i = 0; i + 1 < us.length; i++) {
      const [u0, u1] = [us[i]!, us[i + 1]!];
      for (let j = 0; j < vs.length; j++) {
        const [a, b] = [grid[i]![j]!, grid[i + 1]![j]!];
        const mid = surface.pointAt((u0 + u1) / 2, vs[j]!);
        const quarter = surface.pointAt(u0 + (u1 - u0) / 4, vs[j]!);
        const late = surface.pointAt(u1 - (u1 - u0) / 4, vs[j]!);
        const off = Math.max(chordDeviation(mid, a, b), chordDeviation(quarter, a, mid), chordDeviation(late, mid, b));
        if (off > tolerance || a.distanceTo(b) > maxEdge) {
          splitU.add(i);
          break;
        }
      }
    }
    for (let j = 0; j + 1 < vs.length; j++) {
      const [v0, v1] = [vs[j]!, vs[j + 1]!];
      for (let i = 0; i < us.length; i++) {
        const [a, b] = [grid[i]![j]!, grid[i]![j + 1]!];
        const mid = surface.pointAt(us[i]!, (v0 + v1) / 2);
        const quarter = surface.pointAt(us[i]!, v0 + (v1 - v0) / 4);
        const late = surface.pointAt(us[i]!, v1 - (v1 - v0) / 4);
        const off = Math.max(chordDeviation(mid, a, b), chordDeviation(quarter, a, mid), chordDeviation(late, mid, b));
        if (off > tolerance || a.distanceTo(b) > maxEdge) {
          splitV.add(j);
          break;
        }
      }
    }
    // Twist shows up only across a cell, along the diagonal the triangles use.
    for (let i = 0; i + 1 < us.length; i++) {
      for (let j = 0; j + 1 < vs.length; j++) {
        if (splitU.has(i) && splitV.has(j)) {
          continue;
        }
        const centre = surface.pointAt((us[i]! + us[i + 1]!) / 2, (vs[j]! + vs[j + 1]!) / 2);
        if (chordDeviation(centre, grid[i]![j]!, grid[i + 1]![j + 1]!) > tolerance) {
          splitU.add(i);
          splitV.add(j);
        }
      }
    }
    if (splitU.size === 0 && splitV.size === 0) {
      break;
    }
    if ((us.length + splitU.size) * (vs.length + splitV.size) > maxPoints) {
      break;
    }
    us = refined(us, splitU);
    vs = refined(vs, splitV);
  }
  return { us, vs };
}

export function tessellateSurface(surface: ISurface, options: SurfaceTessellationOptions = {}): SurfaceMesh {
  const { us, vs } = refineGrid(surface, options);
  const points: Vector3[] = [];
  const normals: Vector3[] = [];
  const parameters: SurfaceParameter[] = [];
  for (const u of us) {
    for (const v of vs) {
      points.push(surface.pointAt(u, v));
      normals.push(surface.normalAt(u, v));
      parameters.push({ u, v });
    }
  }
  const stride = vs.length;
  const triangles: number[] = [];
  const addTriangle = (a: number, b: number, c: number) => {
    const [pa, pb, pc] = [points[a]!, points[b]!, points[c]!];
    const edge = Math.max(pa.distanceTo(pb), pb.distanceTo(pc), pc.distanceTo(pa));
    // Triangles that collapse at a pole have no area.
    if (pb.sub(pa).cross(pc.sub(pa)).length > 1e-12 * edge * edge) {
      triangles.push(a, b, c);
    }
  };
  for (let i = 0; i + 1 < us.length; i++) {
    for (let j = 0; j + 1 < vs.length; j++) {
      const a = i * stride + j;
      const b = (i + 1) * stride + j;
      addTriangle(a, b, b + 1);
      addTriangle(a, b + 1, a + 1);
    }
  }
  return { points, normals, parameters, triangles };
}
