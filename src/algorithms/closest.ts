import type { Vector3 } from 'mc-curves';
import { clamp, gridParameters } from '../math/numeric.js';
import type { ClosestSurfacePoint, ISurface } from '../surfaces/ISurface.js';

const GRID_INTERVALS = 32;
const MAX_STARTS = 8;

/** Damped Newton descent of the squared distance, kept inside the domain. */
function refine(surface: ISurface, target: Vector3, u0: number, v0: number): ClosestSurfacePoint {
  const { domainU: du, domainV: dv } = surface;
  const wrapU = surface.isPeriodicU;
  const wrapV = surface.isPeriodicV;
  const limitU = (u: number) => (wrapU ? u : clamp(u, du.start, du.end));
  const limitV = (v: number) => (wrapV ? v : clamp(v, dv.start, dv.end));
  let u = u0;
  let v = v0;
  let p = surface.pointAt(u, v);
  let best = p.sub(target).lengthSquared;
  for (let iter = 0; iter < 40; iter++) {
    const d = p.sub(target);
    const su = surface.derivativeAt(u, v, 1, 0);
    const sv = surface.derivativeAt(u, v, 0, 1);
    const gu = su.dot(d);
    const gv = sv.dot(d);
    const E = su.lengthSquared;
    const F = su.dot(sv);
    const G = sv.lengthSquared;
    let a = E + surface.derivativeAt(u, v, 2, 0).dot(d);
    let b = F + surface.derivativeAt(u, v, 1, 1).dot(d);
    let c = G + surface.derivativeAt(u, v, 0, 2).dot(d);
    if (!(a > 0 && c > 0 && a * c - b * b > 0)) {
      // Not locally convex, so use the always-descending Gauss-Newton matrix.
      [a, b, c] = [E, F, G];
    }
    const damping = 1e-12 * (a + c) + 1e-300;
    a += damping;
    c += damping;
    const det = a * c - b * b;
    let stepU = det > 0 ? -(c * gu - b * gv) / det : 0;
    let stepV = det > 0 ? -(a * gv - b * gu) / det : 0;
    // A variable pressed against its limit stays there, and the other moves alone.
    if (!wrapU && ((u <= du.start && stepU < 0) || (u >= du.end && stepU > 0))) {
      stepU = 0;
      stepV = -gv / c;
    } else if (!wrapV && ((v <= dv.start && stepV < 0) || (v >= dv.end && stepV > 0))) {
      stepV = 0;
      stepU = -gu / a;
    }
    let scale = 1;
    let accepted = false;
    for (let k = 0; k < 12; k++) {
      const nu = limitU(u + scale * stepU);
      const nv = limitV(v + scale * stepV);
      const np = surface.pointAt(nu, nv);
      const f = np.sub(target).lengthSquared;
      if (f < best) {
        const moved = Math.abs(nu - u) / (du.end - du.start) + Math.abs(nv - v) / (dv.end - dv.start);
        u = nu;
        v = nv;
        p = np;
        best = f;
        accepted = moved > 1e-13;
        break;
      }
      scale /= 2;
    }
    if (!accepted) {
      break;
    }
  }
  return { u, v, point: p, distance: Math.sqrt(best) };
}

/** Samples a grid, then refines each local minimum, as there can be several on a curved surface. */
export function closestSurfacePoint(surface: ISurface, target: Vector3): ClosestSurfacePoint {
  const us = gridParameters(surface.breakpointsU(), GRID_INTERVALS);
  const vs = gridParameters(surface.breakpointsV(), GRID_INTERVALS);
  const distance = us.map((u) => vs.map((v) => surface.pointAt(u, v).distanceTo(target)));
  const starts: { i: number; j: number; d: number }[] = [];
  for (let i = 0; i < us.length; i++) {
    for (let j = 0; j < vs.length; j++) {
      const d = distance[i]![j]!;
      let isMinimum = true;
      for (let a = Math.max(0, i - 1); a <= Math.min(us.length - 1, i + 1) && isMinimum; a++) {
        for (let b = Math.max(0, j - 1); b <= Math.min(vs.length - 1, j + 1); b++) {
          if (distance[a]![b]! < d) {
            isMinimum = false;
            break;
          }
        }
      }
      if (isMinimum) {
        starts.push({ i, j, d });
      }
    }
  }
  starts.sort((a, b) => a.d - b.d);
  let best: ClosestSurfacePoint | undefined;
  for (const start of starts.slice(0, MAX_STARTS)) {
    const result = refine(surface, target, us[start.i]!, vs[start.j]!);
    if (!best || result.distance < best.distance) {
      best = result;
    }
  }
  const { u, v, point, distance: d } = best!;
  const span = (t: number, start: number, end: number) => start + ((((t - start) % (end - start)) + (end - start)) % (end - start));
  return {
    u: surface.isPeriodicU ? span(u, surface.domainU.start, surface.domainU.end) : u,
    v: surface.isPeriodicV ? span(v, surface.domainV.start, surface.domainV.end) : v,
    point,
    distance: d,
  };
}
