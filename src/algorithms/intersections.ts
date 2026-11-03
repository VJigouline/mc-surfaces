import { BoundingBox3, Line, Vector3, type ICurve } from 'meercad-curves';
import { bisect, clamp, solveLinear } from '../math/numeric.js';
import type {
  CurveSurfaceIntersection,
  ISurface,
  SurfaceIntersection,
  SurfaceIntersectionOptions,
  SurfaceMesh,
  SurfaceParameter,
  SurfaceSection,
} from '../surfaces/ISurface.js';

const MAX_CURVE_SEEDS = 2000;

/** Axis-aligned box tree over a mesh's triangles. */
class TriangleTree {
  private readonly root: TreeNode;

  constructor(readonly mesh: SurfaceMesh) {
    const count = mesh.triangles.length / 3;
    const items = Array.from({ length: count }, (_, i) => i);
    this.root = this.build(items);
  }

  corner(triangle: number, k: number): number {
    return this.mesh.triangles[triangle * 3 + k]!;
  }

  private boxOf(triangle: number): BoundingBox3 {
    const { points } = this.mesh;
    return BoundingBox3.fromPoints([points[this.corner(triangle, 0)]!, points[this.corner(triangle, 1)]!, points[this.corner(triangle, 2)]!]);
  }

  private build(items: number[]): TreeNode {
    const corners = items.flatMap((t) => [0, 1, 2].map((k) => this.mesh.points[this.corner(t, k)]!));
    const box = corners.length > 0 ? BoundingBox3.fromPoints(corners) : new BoundingBox3(Vector3.ZERO, Vector3.ZERO);
    if (items.length <= 4) {
      return { box, items };
    }
    const size = box.size;
    const axis = size.x >= size.y && size.x >= size.z ? 'x' : size.y >= size.z ? 'y' : 'z';
    const centre = (t: number) => this.boxOf(t).center[axis];
    const sorted = items.map((t) => [centre(t), t] as const).sort((a, b) => a[0] - b[0]);
    const half = sorted.length >> 1;
    return {
      box,
      left: this.build(sorted.slice(0, half).map((e) => e[1])),
      right: this.build(sorted.slice(half).map((e) => e[1])),
    };
  }

  /** Calls `visit` for each triangle whose box may touch `box`. */
  query(box: BoundingBox3, margin: number, visit: (triangle: number) => void): void {
    const walk = (node: TreeNode): void => {
      if (!node.box.intersects(box, margin)) {
        return;
      }
      if (node.items) {
        node.items.forEach(visit);
      } else {
        walk(node.left!);
        walk(node.right!);
      }
    };
    if (this.mesh.triangles.length > 0) {
      walk(this.root);
    }
  }
}

interface TreeNode {
  readonly box: BoundingBox3;
  readonly items?: readonly number[];
  readonly left?: TreeNode;
  readonly right?: TreeNode;
}

interface SegmentHit {
  /** Position along the segment, 0 to 1. */
  readonly s: number;
  /** Barycentric weights of the second and third triangle corners. */
  readonly b1: number;
  readonly b2: number;
}

/** Moller-Trumbore intersection of a segment and a triangle, with a small margin on the edges. */
function segmentTriangle(a: Vector3, b: Vector3, p0: Vector3, p1: Vector3, p2: Vector3, margin = 1e-9): SegmentHit | undefined {
  const dir = b.sub(a);
  const e1 = p1.sub(p0);
  const e2 = p2.sub(p0);
  const h = dir.cross(e2);
  const det = e1.dot(h);
  if (Math.abs(det) <= 1e-14 * dir.length * e1.length * e2.length) {
    return undefined;
  }
  const inv = 1 / det;
  const o = a.sub(p0);
  const b1 = o.dot(h) * inv;
  if (b1 < -margin || b1 > 1 + margin) {
    return undefined;
  }
  const q = o.cross(e1);
  const b2 = dir.dot(q) * inv;
  if (b2 < -margin || b1 + b2 > 1 + margin) {
    return undefined;
  }
  const s = e2.dot(q) * inv;
  if (s < -margin || s > 1 + margin) {
    return undefined;
  }
  return { s: clamp(s, 0, 1), b1, b2 };
}

function meshParameter(mesh: SurfaceMesh, tree: TriangleTree, triangle: number, b1: number, b2: number): SurfaceParameter {
  const p0 = mesh.parameters[tree.corner(triangle, 0)]!;
  const p1 = mesh.parameters[tree.corner(triangle, 1)]!;
  const p2 = mesh.parameters[tree.corner(triangle, 2)]!;
  const b0 = 1 - b1 - b2;
  return { u: b0 * p0.u + b1 * p1.u + b2 * p2.u, v: b0 * p0.v + b1 * p1.v + b2 * p2.v };
}

function limit(surface: ISurface, p: SurfaceParameter): SurfaceParameter {
  const { domainU: du, domainV: dv } = surface;
  return {
    u: surface.isPeriodicU ? p.u : clamp(p.u, du.start, du.end),
    v: surface.isPeriodicV ? p.v : clamp(p.v, dv.start, dv.end),
  };
}

function defaultMeshTolerance(...surfaces: readonly ISurface[]): number {
  return Math.min(...surfaces.map((s) => Math.max(s.boundingBox().diagonal, 1e-300))) / 300;
}

/** Levenberg-Marquardt on `S(u, v) - C(t) = 0`. */
function refineCurveHit(surface: ISurface, curve: ICurve, seed: CurveSurfaceIntersection): CurveSurfaceIntersection {
  const domain = curve.domain;
  let { u, v, t } = seed;
  let p = surface.pointAt(u, v);
  let q = curve.pointAt(t);
  let f = p.sub(q);
  let cost = f.lengthSquared;
  let lambda = 1e-12;
  for (let iter = 0; iter < 40 && cost > 0; iter++) {
    const su = surface.derivativeAt(u, v, 1, 0);
    const sv = surface.derivativeAt(u, v, 0, 1);
    const ct = curve.derivativeAt(t).negate();
    const J = [su, sv, ct];
    const JtJ = J.map((a) => J.map((b) => a.dot(b)));
    const g = J.map((a) => -a.dot(f));
    const trace = JtJ[0]![0]! + JtJ[1]![1]! + JtJ[2]![2]!;
    let improved = false;
    for (let k = 0; k < 8; k++) {
      const m = JtJ.map((row, i) => row.map((value, j) => (i === j ? value + lambda * trace + 1e-300 : value)));
      let step: number[];
      try {
        step = solveLinear(m, g);
      } catch {
        lambda *= 100;
        continue;
      }
      const next = limit(surface, { u: u + step[0]!, v: v + step[1]! });
      const nt = clamp(t + step[2]!, domain.start, domain.end);
      const np = surface.pointAt(next.u, next.v);
      const nq = curve.pointAt(nt);
      const nf = np.sub(nq);
      if (nf.lengthSquared < cost) {
        u = next.u;
        v = next.v;
        t = nt;
        p = np;
        q = nq;
        f = nf;
        cost = nf.lengthSquared;
        lambda = Math.max(lambda / 10, 1e-14);
        improved = true;
        break;
      }
      lambda *= 10;
    }
    if (!improved) {
      break;
    }
  }
  return { u, v, t, point: p, distance: Math.sqrt(cost) };
}

function wrapParameters(surface: ISurface, p: SurfaceParameter): SurfaceParameter {
  const wrap = (x: number, d: { start: number; end: number }) => d.start + ((((x - d.start) % (d.end - d.start)) + (d.end - d.start)) % (d.end - d.start));
  return {
    u: surface.isPeriodicU ? wrap(p.u, surface.domainU) : p.u,
    v: surface.isPeriodicV ? wrap(p.v, surface.domainV) : p.v,
  };
}

export function intersectCurveSurface(surface: ISurface, curve: ICurve, tolerance = 1e-6): CurveSurfaceIntersection[] {
  const box = surface.boundingBox();
  const curveBox = curve.boundingBox();
  const meshTolerance = Math.max(Math.min(box.diagonal, curveBox.diagonal), box.diagonal * 1e-3) / 300;
  if (!box.intersects(curveBox, meshTolerance * 2 + tolerance)) {
    return [];
  }
  const mesh = surface.tessellate({ chordTolerance: meshTolerance });
  const tree = new TriangleTree(mesh);
  const polyline = curve.tessellate({ chordTolerance: meshTolerance });
  const seeds: CurveSurfaceIntersection[] = [];
  const near = 2 * meshTolerance;
  for (let i = 0; i + 1 < polyline.points.length && seeds.length < MAX_CURVE_SEEDS; i++) {
    const [a, b] = [polyline.points[i]!, polyline.points[i + 1]!];
    const [t0, t1] = [polyline.parameters[i]!, polyline.parameters[i + 1]!];
    const mid = a.lerp(b, 0.5);
    tree.query(BoundingBox3.fromPoints([a, b]), near, (triangle) => {
      if (seeds.length >= MAX_CURVE_SEEDS) {
        return;
      }
      const [p0, p1, p2] = [0, 1, 2].map((k) => mesh.points[tree.corner(triangle, k)]!) as [Vector3, Vector3, Vector3];
      const hit = segmentTriangle(a, b, p0, p1, p2, 1e-6);
      if (hit) {
        const uv = meshParameter(mesh, tree, triangle, hit.b1, hit.b2);
        seeds.push({ ...uv, t: t0 + (t1 - t0) * hit.s, point: a, distance: Infinity });
        return;
      }
      // A curve that touches without crossing never pierces the mesh, so also try close passes.
      const normal = p1.sub(p0).cross(p2.sub(p0));
      const area2 = normal.lengthSquared;
      if (area2 === 0) {
        return;
      }
      const height = mid.sub(p0).dot(normal) / Math.sqrt(area2);
      if (Math.abs(height) > near) {
        return;
      }
      const off = mid.sub(p0);
      const b1 = off.cross(p2.sub(p0)).dot(normal) / area2;
      const b2 = p1.sub(p0).cross(off).dot(normal) / area2;
      if (b1 >= 0 && b2 >= 0 && b1 + b2 <= 1) {
        const uv = meshParameter(mesh, tree, triangle, b1, b2);
        seeds.push({ ...uv, t: (t0 + t1) / 2, point: mid, distance: Infinity });
      }
    });
  }
  const merge = Math.max(10 * tolerance, 1e-4 * Math.min(box.diagonal, curveBox.diagonal));
  const results: CurveSurfaceIntersection[] = [];
  for (const seed of seeds) {
    const hit = refineCurveHit(surface, curve, seed);
    if (hit.distance > tolerance) {
      continue;
    }
    const twin = results.findIndex((r) => r.point.distanceTo(hit.point) <= merge);
    if (twin < 0) {
      results.push(hit);
    } else if (hit.distance < results[twin]!.distance) {
      results[twin] = hit;
    }
  }
  return results
    .map((r) => ({ ...r, ...wrapParameters(surface, r) }))
    .sort((a, b) => a.t - b.t);
}

export function intersectLineSurface(surface: ISurface, origin: Vector3, direction: Vector3, tolerance = 1e-6): CurveSurfaceIntersection[] {
  const dir = direction.normalized();
  if (dir.lengthSquared === 0) {
    throw new RangeError('Line direction must not be zero.');
  }
  const box = surface.boundingBox();
  const centre = origin.addScaled(dir, box.center.sub(origin).dot(dir));
  const reach = box.diagonal + centre.distanceTo(box.center) + 1e-9;
  const line = new Line(centre.addScaled(dir, -reach), centre.addScaled(dir, reach));
  return intersectCurveSurface(surface, line, tolerance).map((hit) => ({ ...hit, t: hit.point.sub(origin).dot(dir) }));
}

interface Endpoint<T> {
  readonly point: Vector3;
  readonly data: T;
}

interface Chain<T> {
  readonly nodes: Endpoint<T>[];
  readonly closed: boolean;
}

/** Joins segments that share end points into polylines. Branch points end a polyline. */
function chainSegments<T>(segments: readonly [Endpoint<T>, Endpoint<T>][], tolerance: number): Chain<T>[] {
  const nodes: Endpoint<T>[] = [];
  const cells = new Map<string, number[]>();
  const nodeOf = (end: Endpoint<T>): number => {
    const key = [end.point.x, end.point.y, end.point.z].map((c) => Math.floor(c / tolerance));
    for (let a = -1; a <= 1; a++) {
      for (let b = -1; b <= 1; b++) {
        for (let c = -1; c <= 1; c++) {
          for (const id of cells.get(`${key[0]! + a},${key[1]! + b},${key[2]! + c}`) ?? []) {
            if (nodes[id]!.point.distanceTo(end.point) <= tolerance) {
              return id;
            }
          }
        }
      }
    }
    const id = nodes.length;
    nodes.push(end);
    const name = key.join(',');
    cells.set(name, [...(cells.get(name) ?? []), id]);
    return id;
  };
  const links = new Map<number, number[]>();
  const seen = new Set<string>();
  for (const [a, b] of segments) {
    const [i, j] = [nodeOf(a), nodeOf(b)];
    const name = i < j ? `${i}/${j}` : `${j}/${i}`;
    if (i === j || seen.has(name)) {
      continue;
    }
    seen.add(name);
    links.set(i, [...(links.get(i) ?? []), j]);
    links.set(j, [...(links.get(j) ?? []), i]);
  }
  const used = new Set<string>();
  const edge = (i: number, j: number) => (i < j ? `${i}/${j}` : `${j}/${i}`);
  const chains: Chain<T>[] = [];
  const walk = (start: number): void => {
    for (const first of links.get(start)!) {
      if (used.has(edge(start, first))) {
        continue;
      }
      used.add(edge(start, first));
      const ids = [start, first];
      let current = first;
      while (links.get(current)!.length === 2) {
        const next = links.get(current)!.find((n) => !used.has(edge(current, n)));
        if (next === undefined) {
          break;
        }
        used.add(edge(current, next));
        ids.push(next);
        current = next;
      }
      const closed = current === start;
      if (closed) {
        ids.pop();
      }
      chains.push({ nodes: ids.map((id) => nodes[id]!), closed });
    }
  };
  // Open ends and branch points first, then whatever is left must be a loop.
  for (const [id, neighbours] of links) {
    if (neighbours.length !== 2) {
      walk(id);
    }
  }
  for (const id of links.keys()) {
    walk(id);
  }
  return chains;
}

export function intersectPlaneSurface(surface: ISurface, origin: Vector3, normal: Vector3, options: SurfaceIntersectionOptions = {}): SurfaceSection[] {
  const n = normal.normalized();
  if (n.lengthSquared === 0) {
    throw new RangeError('Plane normal must not be zero.');
  }
  const meshTolerance = options.meshTolerance ?? defaultMeshTolerance(surface);
  const mesh = surface.tessellate({ chordTolerance: meshTolerance });
  // Heights within rounding noise of zero count as positive, on every side of a seam alike. Otherwise
  // sin(2 pi) being slightly negative while sin(0) is zero would lose a section running along the seam.
  const noise = 1e-13 * Math.max(surface.boundingBox().diagonal, origin.length);
  const height = (p: Vector3) => {
    const h = p.sub(origin).dot(n);
    return Math.abs(h) <= noise ? noise : h;
  };
  const heights = mesh.points.map(height);
  const crossings = new Map<string, Endpoint<SurfaceParameter>>();
  const crossing = (a: number, b: number): Endpoint<SurfaceParameter> => {
    const key = a < b ? `${a}/${b}` : `${b}/${a}`;
    let result = crossings.get(key);
    if (!result) {
      const [pa, pb] = [mesh.parameters[a]!, mesh.parameters[b]!];
      const at = (s: number) => ({ u: pa.u + (pb.u - pa.u) * s, v: pa.v + (pb.v - pa.v) * s });
      const s = bisect((x) => {
        const q = at(x);
        return height(surface.pointAt(q.u, q.v));
      }, 0, 1);
      const q = at(s);
      result = { point: surface.pointAt(q.u, q.v), data: q };
      crossings.set(key, result);
    }
    return result;
  };
  const segments: [Endpoint<SurfaceParameter>, Endpoint<SurfaceParameter>][] = [];
  for (let t = 0; t < mesh.triangles.length; t += 3) {
    const corners = [mesh.triangles[t]!, mesh.triangles[t + 1]!, mesh.triangles[t + 2]!];
    const ends: Endpoint<SurfaceParameter>[] = [];
    for (let k = 0; k < 3; k++) {
      const [a, b] = [corners[k]!, corners[(k + 1) % 3]!];
      if (heights[a]! >= 0 !== heights[b]! >= 0) {
        ends.push(crossing(a, b));
      }
    }
    if (ends.length === 2) {
      segments.push([ends[0]!, ends[1]!]);
    }
  }
  const merge = 1e-8 * Math.max(surface.boundingBox().diagonal, 1e-300);
  return chainSegments(segments, merge).map((chain) => ({
    points: chain.nodes.map((node) => node.point),
    parameters: chain.nodes.map((node) => node.data),
    closed: chain.closed,
  }));
}

interface PairParameter {
  readonly a: SurfaceParameter;
  readonly b: SurfaceParameter;
}

/** Minimum-norm Gauss-Newton on `A(u, v) - B(s, t) = 0`, which moves both parameters as little as it can. */
function refinePair(a: ISurface, b: ISurface, start: PairParameter): { pair: PairParameter; point: Vector3; gap: number } {
  let pa = start.a;
  let pb = start.b;
  let f = a.pointAt(pa.u, pa.v).sub(b.pointAt(pb.u, pb.v));
  let cost = f.lengthSquared;
  for (let iter = 0; iter < 12 && cost > 0; iter++) {
    const J = [a.derivativeAt(pa.u, pa.v, 1, 0), a.derivativeAt(pa.u, pa.v, 0, 1), b.derivativeAt(pb.u, pb.v, 1, 0).negate(), b.derivativeAt(pb.u, pb.v, 0, 1).negate()];
    const row = (axis: 'x' | 'y' | 'z') => J.map((c) => c[axis]);
    const rows = [row('x'), row('y'), row('z')];
    const JJt = rows.map((r) => rows.map((s) => r.reduce((sum, value, i) => sum + value * s[i]!, 0)));
    const trace = JJt[0]![0]! + JJt[1]![1]! + JJt[2]![2]!;
    const damped = JJt.map((r, i) => r.map((value, j) => (i === j ? value + 1e-10 * trace + 1e-300 : value)));
    let y: number[];
    try {
      y = solveLinear(damped, [-f.x, -f.y, -f.z]);
    } catch {
      break;
    }
    const step = [0, 1, 2, 3].map((c) => rows[0]![c]! * y[0]! + rows[1]![c]! * y[1]! + rows[2]![c]! * y[2]!);
    const na = limit(a, { u: pa.u + step[0]!, v: pa.v + step[1]! });
    const nb = limit(b, { u: pb.u + step[2]!, v: pb.v + step[3]! });
    const nf = a.pointAt(na.u, na.v).sub(b.pointAt(nb.u, nb.v));
    if (!(nf.lengthSquared < cost)) {
      break;
    }
    pa = na;
    pb = nb;
    f = nf;
    cost = nf.lengthSquared;
  }
  const pointA = a.pointAt(pa.u, pa.v);
  return {
    pair: { a: wrapParameters(a, pa), b: wrapParameters(b, pb) },
    point: pointA.lerp(b.pointAt(pb.u, pb.v), 0.5),
    gap: Math.sqrt(cost),
  };
}

/** Triangle pairs of two meshes give segments, which are chained and then pulled onto both surfaces. */
export function intersectSurfaces(a: ISurface, b: ISurface, options: SurfaceIntersectionOptions = {}): SurfaceIntersection[] {
  const tolerance = options.tolerance ?? 1e-6;
  const boxA = a.boundingBox();
  const boxB = b.boundingBox();
  if (!boxA.intersects(boxB, tolerance)) {
    return [];
  }
  const meshTolerance = options.meshTolerance ?? defaultMeshTolerance(a, b);
  const meshA = a.tessellate({ chordTolerance: meshTolerance });
  const meshB = b.tessellate({ chordTolerance: meshTolerance });
  const treeA = new TriangleTree(meshA);
  const treeB = new TriangleTree(meshB);
  const merge = 1e-8 * Math.max(boxA.diagonal, boxB.diagonal);
  const segments: [Endpoint<PairParameter>, Endpoint<PairParameter>][] = [];
  const corners = (mesh: SurfaceMesh, tree: TriangleTree, t: number) =>
    [0, 1, 2].map((k) => mesh.points[tree.corner(t, k)]!) as [Vector3, Vector3, Vector3];
  for (let ta = 0; ta < meshA.triangles.length / 3; ta++) {
    const triA = corners(meshA, treeA, ta);
    treeB.query(BoundingBox3.fromPoints(triA), merge, (tb) => {
      const triB = corners(meshB, treeB, tb);
      const ends: Endpoint<PairParameter>[] = [];
      for (let k = 0; k < 3; k++) {
        const [i, j] = [k, (k + 1) % 3];
        const hitA = segmentTriangle(triA[i]!, triA[j]!, ...triB);
        if (hitA) {
          const pi = meshA.parameters[treeA.corner(ta, i)]!;
          const pj = meshA.parameters[treeA.corner(ta, j)]!;
          ends.push({
            point: triA[i]!.lerp(triA[j]!, hitA.s),
            data: {
              a: { u: pi.u + (pj.u - pi.u) * hitA.s, v: pi.v + (pj.v - pi.v) * hitA.s },
              b: meshParameter(meshB, treeB, tb, hitA.b1, hitA.b2),
            },
          });
        }
        const hitB = segmentTriangle(triB[i]!, triB[j]!, ...triA);
        if (hitB) {
          const pi = meshB.parameters[treeB.corner(tb, i)]!;
          const pj = meshB.parameters[treeB.corner(tb, j)]!;
          ends.push({
            point: triB[i]!.lerp(triB[j]!, hitB.s),
            data: {
              a: meshParameter(meshA, treeA, ta, hitB.b1, hitB.b2),
              b: { u: pi.u + (pj.u - pi.u) * hitB.s, v: pi.v + (pj.v - pi.v) * hitB.s },
            },
          });
        }
      }
      if (ends.length >= 2) {
        // Near a shared edge the same crossing can appear twice, so take the two furthest apart.
        let [first, second, widest] = [ends[0]!, ends[1]!, -1];
        for (let i = 0; i < ends.length; i++) {
          for (let j = i + 1; j < ends.length; j++) {
            const d = ends[i]!.point.distanceTo(ends[j]!.point);
            if (d > widest) {
              [first, second, widest] = [ends[i]!, ends[j]!, d];
            }
          }
        }
        segments.push([first, second]);
      }
    });
  }
  const results: SurfaceIntersection[] = [];
  for (const chain of chainSegments(segments, merge)) {
    const points: Vector3[] = [];
    const parameters: SurfaceParameter[] = [];
    const parametersOther: SurfaceParameter[] = [];
    for (const node of chain.nodes) {
      const refined = refinePair(a, b, node.data);
      if (refined.gap > tolerance) {
        continue;
      }
      if (points.length > 0 && points[points.length - 1]!.distanceTo(refined.point) <= merge) {
        continue;
      }
      points.push(refined.point);
      parameters.push(refined.pair.a);
      parametersOther.push(refined.pair.b);
    }
    if (points.length >= 2) {
      results.push({ points, parameters, parametersOther, closed: chain.closed });
    }
  }
  return results;
}
