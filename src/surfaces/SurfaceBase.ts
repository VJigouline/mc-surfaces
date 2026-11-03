import { BoundingBox3, CurveBase, Transform3, Vector3, type CurveDomain, type CurveJSON, type ICurve } from 'meercad-curves';
import { closestSurfacePoint } from '../algorithms/closest.js';
import { intersectCurveSurface, intersectLineSurface, intersectPlaneSurface, intersectSurfaces } from '../algorithms/intersections.js';
import { tessellateSurface } from '../algorithms/tessellate.js';
import { clamp, finiteDifferenceWeights, GEOMETRIC_TOLERANCE, GL_NODES, GL_WEIGHTS, goldenMaximum, gridParameters } from '../math/numeric.js';
import type {
  ClosestSurfacePoint,
  CurveSurfaceIntersection,
  FirstFundamentalForm,
  ISurface,
  PrincipalCurvatures,
  SecondFundamentalForm,
  SurfaceDomain,
  SurfaceFrame,
  SurfaceIntersection,
  SurfaceIntersectionOptions,
  SurfaceJSON,
  SurfaceMesh,
  SurfaceSection,
  SurfaceTessellationOptions,
} from './ISurface.js';
import type { NurbsSurface } from './NurbsSurface.js';

type NurbsApproximator = (surface: ISurface, tolerance: number) => NurbsSurface;

let nurbsApproximator: NurbsApproximator | undefined;

/**
 * Lets `NurbsSurface` give the base class a way to build it. A direct import would be a circular
 * module dependency, since `NurbsSurface` extends `SurfaceBase`.
 */
export function registerNurbsApproximator(approximator: NurbsApproximator): void {
  nurbsApproximator = approximator;
}

/** Below this fraction of the typical speed, a parameter derivative counts as vanished. */
const SINGULAR_SPEED = 1e-9;
/** How far inside the domain to measure curvature at a singular point, as a fraction of the span. */
const SINGULAR_NUDGE = 1e-4;
const AREA_MAX_DEPTH = 6;

function gaussCell(f: (u: number, v: number) => number, u0: number, u1: number, v0: number, v1: number): number {
  const mu = (u0 + u1) / 2;
  const hu = (u1 - u0) / 2;
  const mv = (v0 + v1) / 2;
  const hv = (v1 - v0) / 2;
  let sum = 0;
  for (let i = 0; i < GL_NODES.length; i++) {
    const du = hu * GL_NODES[i]!;
    for (let j = 0; j < GL_NODES.length; j++) {
      const dv = hv * GL_NODES[j]!;
      sum += GL_WEIGHTS[i]! * GL_WEIGHTS[j]! * (f(mu - du, mv - dv) + f(mu - du, mv + dv) + f(mu + du, mv - dv) + f(mu + du, mv + dv));
    }
  }
  return sum * hu * hv;
}

/** Adaptive integration over a rectangle, comparing each cell with its four quarters. */
function integrateCell(f: (u: number, v: number) => number, u0: number, u1: number, v0: number, v1: number, whole: number, floor: number, depth: number): number {
  const um = (u0 + u1) / 2;
  const vm = (v0 + v1) / 2;
  const a = gaussCell(f, u0, um, v0, vm);
  const b = gaussCell(f, um, u1, v0, vm);
  const c = gaussCell(f, u0, um, vm, v1);
  const d = gaussCell(f, um, u1, vm, v1);
  const sum = a + b + c + d;
  if (depth <= 0 || Math.abs(sum - whole) <= floor) {
    return sum;
  }
  return (
    integrateCell(f, u0, um, v0, vm, a, floor / 2, depth - 1) +
    integrateCell(f, um, u1, v0, vm, b, floor / 2, depth - 1) +
    integrateCell(f, u0, um, vm, v1, c, floor / 2, depth - 1) +
    integrateCell(f, um, u1, vm, v1, d, floor / 2, depth - 1)
  );
}

/** The breakpoints inside [lo, hi], plus both limits. */
function clipBreaks(breaks: readonly number[], lo: number, hi: number): number[] {
  return [lo, ...breaks.filter((b) => b > lo + 1e-12 && b < hi - 1e-12), hi];
}

/**
 * Numeric implementation of everything in `ISurface`, given only `evaluate`.
 * Subclasses override methods where a closed form is faster or more accurate.
 */
export abstract class SurfaceBase implements ISurface {
  abstract readonly type: string;
  abstract readonly domainU: SurfaceDomain;
  abstract readonly domainV: SurfaceDomain;
  private cachedArea: number | undefined;
  private cachedBox: BoundingBox3 | undefined;
  private cachedSpeeds: [number, number] | undefined;

  /** Point at parameters that are already inside the domain. */
  protected abstract evaluate(u: number, v: number): Vector3;
  abstract toJSON(): SurfaceJSON;

  /** Derivative at parameters already inside the domain. The default is numeric, up to order 4 per parameter. */
  protected evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.numericDerivative(u, v, orderU, orderV);
  }

  get degreeU(): number | undefined {
    return undefined;
  }

  get degreeV(): number | undefined {
    return undefined;
  }

  get isPeriodicU(): boolean {
    return false;
  }

  get isPeriodicV(): boolean {
    return false;
  }

  private edgeTolerance(p: Vector3): number {
    return GEOMETRIC_TOLERANCE + 1e-12 * p.length;
  }

  private edgesCoincide(a: (s: number) => Vector3, b: (s: number) => Vector3): boolean {
    for (let i = 0; i <= 8; i++) {
      const p = a(i / 8);
      if (p.distanceTo(b(i / 8)) > this.edgeTolerance(p)) {
        return false;
      }
    }
    return true;
  }

  private atU(s: number): number {
    return this.domainU.start + (this.domainU.end - this.domainU.start) * s;
  }

  private atV(s: number): number {
    return this.domainV.start + (this.domainV.end - this.domainV.start) * s;
  }

  get isClosedU(): boolean {
    const { start, end } = this.domainU;
    return this.edgesCoincide((s) => this.evaluate(start, this.atV(s)), (s) => this.evaluate(end, this.atV(s)));
  }

  get isClosedV(): boolean {
    const { start, end } = this.domainV;
    return this.edgesCoincide((s) => this.evaluate(this.atU(s), start), (s) => this.evaluate(this.atU(s), end));
  }

  get enclosesVolume(): boolean {
    const { domainU: du, domainV: dv } = this;
    const collapsed = (edge: (s: number) => Vector3) => this.edgesCoincide(edge, () => edge(0));
    const uOk = this.isClosedU || (collapsed((s) => this.evaluate(du.start, this.atV(s))) && collapsed((s) => this.evaluate(du.end, this.atV(s))));
    const vOk = this.isClosedV || (collapsed((s) => this.evaluate(this.atU(s), dv.start)) && collapsed((s) => this.evaluate(this.atU(s), dv.end)));
    return uOk && vOk;
  }

  get area(): number {
    this.cachedArea ??= this.areaBetween(this.domainU.start, this.domainU.end, this.domainV.start, this.domainV.end);
    return this.cachedArea;
  }

  /** Clamps a parameter to its domain, or wraps it in a periodic direction. */
  private normalize(t: number, domain: SurfaceDomain, periodic: boolean): number {
    if (Number.isNaN(t)) {
      throw new RangeError('Surface parameter must not be NaN.');
    }
    const { start, end } = domain;
    if (t >= start && t <= end) {
      return t;
    }
    if (periodic) {
      const span = end - start;
      return start + ((((t - start) % span) + span) % span);
    }
    return clamp(t, start, end);
  }

  protected normalizeU(u: number): number {
    return this.normalize(u, this.domainU, this.isPeriodicU);
  }

  protected normalizeV(v: number): number {
    return this.normalize(v, this.domainV, this.isPeriodicV);
  }

  pointAt(u: number, v: number): Vector3 {
    return this.evaluate(this.normalizeU(u), this.normalizeV(v));
  }

  derivativeAt(u: number, v: number, orderU = 1, orderV = 0): Vector3 {
    if (!Number.isInteger(orderU) || !Number.isInteger(orderV) || orderU < 0 || orderV < 0) {
      throw new RangeError('Derivative orders must be non-negative integers.');
    }
    const nu = this.normalizeU(u);
    const nv = this.normalizeV(v);
    return orderU === 0 && orderV === 0 ? this.evaluate(nu, nv) : this.evaluateDerivative(nu, nv, orderU, orderV);
  }

  /** The seven nodes for a derivative at `t`, shifted to stay inside its smooth span. */
  private stencil(t: number, breaks: readonly number[], order: number): { nodes: number[]; weights: number[] } {
    if (order === 0) {
      return { nodes: [t], weights: [1] };
    }
    let lo = breaks[0]!;
    let hi = breaks[breaks.length - 1]!;
    for (let i = 0; i + 1 < breaks.length; i++) {
      if (t >= breaks[i]! && t < breaks[i + 1]!) {
        lo = breaks[i]!;
        hi = breaks[i + 1]!;
        break;
      }
    }
    if (t >= hi && breaks.length >= 2) {
      lo = breaks[breaks.length - 2]!;
    }
    const h = (hi - lo) * 1e-3;
    const first = clamp(t - 3 * h, lo, hi - 6 * h);
    const nodes = Array.from({ length: 7 }, (_, i) => first + i * h);
    return { nodes, weights: finiteDifferenceWeights(t, nodes, order) };
  }

  protected numericDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    if (orderU > 4 || orderV > 4) {
      throw new RangeError(`${this.type} supports derivatives up to order 4 per parameter.`);
    }
    const su = this.stencil(u, this.breakpointsU(), orderU);
    const sv = this.stencil(v, this.breakpointsV(), orderV);
    let sum = Vector3.ZERO;
    su.nodes.forEach((nu, i) => {
      sv.nodes.forEach((nv, j) => {
        sum = sum.addScaled(this.evaluate(nu, nv), su.weights[i]! * sv.weights[j]!);
      });
    });
    return sum;
  }

  /** Typical sizes of `Su` and `Sv`, the scale for deciding that one has vanished. */
  private typicalSpeeds(): [number, number] {
    if (!this.cachedSpeeds) {
      let [top, side] = [0, 0];
      for (let i = 0; i < 5; i++) {
        for (let j = 0; j < 5; j++) {
          const u = this.atU((i + 0.5) / 5);
          const v = this.atV((j + 0.5) / 5);
          top = Math.max(top, this.derivativeAt(u, v, 1, 0).length);
          side = Math.max(side, this.derivativeAt(u, v, 0, 1).length);
        }
      }
      this.cachedSpeeds = [top, side];
    }
    return this.cachedSpeeds;
  }

  private isSingular(su: Vector3, sv: Vector3): boolean {
    const [typicalU, typicalV] = this.typicalSpeeds();
    if (su.length <= SINGULAR_SPEED * typicalU || sv.length <= SINGULAR_SPEED * typicalV) {
      return true;
    }
    return su.cross(sv).length <= 1e-10 * su.length * sv.length;
  }

  /** Moves a singular point slightly into the domain. */
  private nudged(u: number, v: number): [number, number] {
    const inward = (t: number, d: SurfaceDomain) => t + (d.end - d.start) * SINGULAR_NUDGE * (t > (d.start + d.end) / 2 ? -1 : 1);
    const su = this.derivativeAt(u, v, 1, 0);
    const sv = this.derivativeAt(u, v, 0, 1);
    const [typicalU, typicalV] = this.typicalSpeeds();
    const flatU = su.length <= SINGULAR_SPEED * typicalU;
    const flatV = sv.length <= SINGULAR_SPEED * typicalV;
    if (flatU && !flatV) {
      return [u, inward(v, this.domainV)];
    }
    if (flatV && !flatU) {
      return [inward(u, this.domainU), v];
    }
    return [inward(u, this.domainU), inward(v, this.domainV)];
  }

  tangentUAt(u: number, v: number): Vector3 {
    const su = this.derivativeAt(u, v, 1, 0);
    if (su.length > SINGULAR_SPEED * this.typicalSpeeds()[0]) {
      return su.normalized();
    }
    const [nu, nv] = this.nudged(u, v);
    return this.derivativeAt(nu, nv, 1, 0).normalized();
  }

  tangentVAt(u: number, v: number): Vector3 {
    const sv = this.derivativeAt(u, v, 0, 1);
    if (sv.length > SINGULAR_SPEED * this.typicalSpeeds()[1]) {
      return sv.normalized();
    }
    const [nu, nv] = this.nudged(u, v);
    return this.derivativeAt(nu, nv, 0, 1).normalized();
  }

  normalAt(u: number, v: number): Vector3 {
    const su = this.derivativeAt(u, v, 1, 0);
    const sv = this.derivativeAt(u, v, 0, 1);
    if (!this.isSingular(su, sv)) {
      return su.cross(sv).normalized();
    }
    // Near an edge where Su vanishes, Su is about (v - v0) Suv, so the limit is exact. Likewise for Sv.
    const nv = this.normalizeV(v);
    const nu = this.normalizeU(u);
    const [typicalU, typicalV] = this.typicalSpeeds();
    let limit = Vector3.ZERO;
    if (su.length <= SINGULAR_SPEED * typicalU) {
      const sign = nv > (this.domainV.start + this.domainV.end) / 2 ? -1 : 1;
      limit = this.derivativeAt(nu, nv, 1, 1).cross(sv).scale(sign);
    } else if (sv.length <= SINGULAR_SPEED * typicalV) {
      const sign = nu > (this.domainU.start + this.domainU.end) / 2 ? -1 : 1;
      limit = su.cross(this.derivativeAt(nu, nv, 1, 1)).scale(sign);
    }
    if (limit.length > 1e-12 * Math.max(typicalU, typicalV) ** 2) {
      return limit.normalized();
    }
    const [mu, mv] = this.nudged(nu, nv);
    return this.derivativeAt(mu, mv, 1, 0).cross(this.derivativeAt(mu, mv, 0, 1)).normalized();
  }

  frameAt(u: number, v: number): SurfaceFrame {
    return {
      origin: this.pointAt(u, v),
      tangentU: this.tangentUAt(u, v),
      tangentV: this.tangentVAt(u, v),
      normal: this.normalAt(u, v),
    };
  }

  firstFundamentalFormAt(u: number, v: number): FirstFundamentalForm {
    const su = this.derivativeAt(u, v, 1, 0);
    const sv = this.derivativeAt(u, v, 0, 1);
    return { E: su.lengthSquared, F: su.dot(sv), G: sv.lengthSquared };
  }

  secondFundamentalFormAt(u: number, v: number): SecondFundamentalForm {
    const n = this.normalAt(u, v);
    return {
      L: this.derivativeAt(u, v, 2, 0).dot(n),
      M: this.derivativeAt(u, v, 1, 1).dot(n),
      N: this.derivativeAt(u, v, 0, 2).dot(n),
    };
  }

  principalCurvaturesAt(u: number, v: number): PrincipalCurvatures {
    let [au, av] = [this.normalizeU(u), this.normalizeV(v)];
    if (this.isSingular(this.derivativeAt(au, av, 1, 0), this.derivativeAt(au, av, 0, 1))) {
      [au, av] = this.nudged(au, av);
    }
    const su = this.derivativeAt(au, av, 1, 0);
    const sv = this.derivativeAt(au, av, 0, 1);
    const n = this.normalAt(au, av);
    const { E, F, G } = this.firstFundamentalFormAt(au, av);
    const { L, M, N } = this.secondFundamentalFormAt(au, av);
    const det = E * G - F * F;
    const K = (L * N - M * M) / det;
    const H = (E * N - 2 * F * M + G * L) / (2 * det);
    const root = Math.sqrt(Math.max(H * H - K, 0));
    const k1 = H + root;
    const k2 = H - root;
    // An eigenvector of the shape operator, from whichever equation is better conditioned.
    const a1 = L - k1 * E;
    const b1 = M - k1 * F;
    const a2 = M - k1 * F;
    const b2 = N - k1 * G;
    const direction =
      Math.hypot(a1, b1) >= Math.hypot(a2, b2) ? su.scale(-b1).addScaled(sv, a1) : su.scale(-b2).addScaled(sv, a2);
    const direction1 = direction.length > 1e-9 * (su.length * sv.length * (Math.abs(k1) + 1)) ? direction.normalized() : su.normalized();
    return { k1, k2, direction1, direction2: n.cross(direction1).normalized() };
  }

  gaussianCurvatureAt(u: number, v: number): number {
    const { k1, k2 } = this.principalCurvaturesAt(u, v);
    return k1 * k2;
  }

  meanCurvatureAt(u: number, v: number): number {
    const { k1, k2 } = this.principalCurvaturesAt(u, v);
    return (k1 + k2) / 2;
  }

  breakpointsU(): number[] {
    return [this.domainU.start, this.domainU.end];
  }

  breakpointsV(): number[] {
    return [this.domainV.start, this.domainV.end];
  }

  isoCurveU(u: number): ICurve {
    return new IsoCurve(this, 'u', this.normalizeU(u));
  }

  isoCurveV(v: number): ICurve {
    return new IsoCurve(this, 'v', this.normalizeV(v));
  }

  boundaryCurves(): [ICurve, ICurve, ICurve, ICurve] {
    return [
      this.isoCurveV(this.domainV.start),
      this.isoCurveU(this.domainU.end),
      this.isoCurveV(this.domainV.end),
      this.isoCurveU(this.domainU.start),
    ];
  }

  /** Integrates over the parameter rectangle, one smooth cell at a time. */
  private integrate(f: (u: number, v: number) => number, u0: number, u1: number, v0: number, v1: number): number {
    const us = clipBreaks(this.breakpointsU(), u0, u1);
    const vs = clipBreaks(this.breakpointsV(), v0, v1);
    let coarse = 0;
    const cells: [number, number, number, number, number][] = [];
    for (let i = 0; i + 1 < us.length; i++) {
      for (let j = 0; j + 1 < vs.length; j++) {
        const whole = gaussCell(f, us[i]!, us[i + 1]!, vs[j]!, vs[j + 1]!);
        coarse += Math.abs(whole);
        cells.push([us[i]!, us[i + 1]!, vs[j]!, vs[j + 1]!, whole]);
      }
    }
    // Relative, so that the work does not depend on the model's size.
    const floor = 1e-10 * coarse;
    return cells.reduce((sum, [a, b, c, d, whole]) => sum + integrateCell(f, a, b, c, d, whole, floor, AREA_MAX_DEPTH), 0);
  }

  areaBetween(u0: number, u1: number, v0: number, v1: number): number {
    const [ua, ub] = [clamp(Math.min(u0, u1), this.domainU.start, this.domainU.end), clamp(Math.max(u0, u1), this.domainU.start, this.domainU.end)];
    const [va, vb] = [clamp(Math.min(v0, v1), this.domainV.start, this.domainV.end), clamp(Math.max(v0, v1), this.domainV.start, this.domainV.end)];
    if (!(ub > ua) || !(vb > va)) {
      return 0;
    }
    return this.integrate((u, v) => this.evaluateDerivative(u, v, 1, 0).cross(this.evaluateDerivative(u, v, 0, 1)).length, ua, ub, va, vb);
  }

  enclosedVolume(): number | undefined {
    if (!this.enclosesVolume) {
      return undefined;
    }
    const c = this.boundingBox().center;
    const { domainU: du, domainV: dv } = this;
    const flux = this.integrate(
      (u, v) => this.evaluate(u, v).sub(c).dot(this.evaluateDerivative(u, v, 1, 0).cross(this.evaluateDerivative(u, v, 0, 1))),
      du.start, du.end, dv.start, dv.end,
    );
    return Math.abs(flux) / 3;
  }

  centroid(): Vector3 {
    const { domainU: du, domainV: dv } = this;
    const moment = (axis: 'x' | 'y' | 'z') =>
      this.integrate(
        (u, v) => this.evaluate(u, v)[axis] * this.evaluateDerivative(u, v, 1, 0).cross(this.evaluateDerivative(u, v, 0, 1)).length,
        du.start, du.end, dv.start, dv.end,
      );
    return new Vector3(moment('x'), moment('y'), moment('z')).scale(1 / this.area);
  }

  sample(countU: number, countV: number): Vector3[][] {
    if (!Number.isInteger(countU) || !Number.isInteger(countV) || countU < 2 || countV < 2) {
      throw new RangeError('Sample counts must be integers of at least 2.');
    }
    return Array.from({ length: countU }, (_, i) =>
      Array.from({ length: countV }, (_, j) => this.evaluate(this.atU(i / (countU - 1)), this.atV(j / (countV - 1)))),
    );
  }

  tessellate(options?: SurfaceTessellationOptions): SurfaceMesh {
    return tessellateSurface(this, options);
  }

  boundingBox(): BoundingBox3 {
    if (!this.cachedBox) {
      const us = gridParameters(this.breakpointsU(), 24);
      const vs = gridParameters(this.breakpointsV(), 24);
      const grid = us.map((u) => vs.map((v) => this.evaluate(u, v)));
      const extremes: Vector3[] = grid.flat();
      for (const axis of ['x', 'y', 'z'] as const) {
        for (const sign of [1, -1]) {
          const peaks: { i: number; j: number; value: number }[] = [];
          for (let i = 0; i < us.length; i++) {
            for (let j = 0; j < vs.length; j++) {
              const value = sign * grid[i]![j]![axis];
              let peak = true;
              for (let a = Math.max(0, i - 1); a <= Math.min(us.length - 1, i + 1) && peak; a++) {
                for (let b = Math.max(0, j - 1); b <= Math.min(vs.length - 1, j + 1); b++) {
                  if (sign * grid[a]![b]![axis] > value) {
                    peak = false;
                    break;
                  }
                }
              }
              if (peak) {
                peaks.push({ i, j, value });
              }
            }
          }
          peaks.sort((a, b) => b.value - a.value);
          for (const { i, j } of peaks.slice(0, 3)) {
            let [u, v] = [us[i]!, vs[j]!];
            const [uLo, uHi] = [us[Math.max(i - 1, 0)]!, us[Math.min(i + 1, us.length - 1)]!];
            const [vLo, vHi] = [vs[Math.max(j - 1, 0)]!, vs[Math.min(j + 1, vs.length - 1)]!];
            for (let round = 0; round < 3; round++) {
              u = goldenMaximum((x) => sign * this.evaluate(x, v)[axis], uLo, uHi);
              v = goldenMaximum((x) => sign * this.evaluate(u, x)[axis], vLo, vHi);
            }
            extremes.push(this.evaluate(u, v));
          }
        }
      }
      this.cachedBox = BoundingBox3.fromPoints(extremes);
    }
    return this.cachedBox;
  }

  closestPoint(point: Vector3): ClosestSurfacePoint {
    return closestSurfacePoint(this, point);
  }

  distanceTo(point: Vector3): number {
    return this.closestPoint(point).distance;
  }

  intersectCurve(curve: ICurve, tolerance?: number): CurveSurfaceIntersection[] {
    return intersectCurveSurface(this, curve, tolerance);
  }

  intersectLine(origin: Vector3, direction: Vector3, tolerance?: number): CurveSurfaceIntersection[] {
    return intersectLineSurface(this, origin, direction, tolerance);
  }

  intersectPlane(origin: Vector3, normal: Vector3, options?: SurfaceIntersectionOptions): SurfaceSection[] {
    return intersectPlaneSurface(this, origin, normal, options);
  }

  intersect(other: ISurface, options?: SurfaceIntersectionOptions): SurfaceIntersection[] {
    return intersectSurfaces(this, other, options);
  }

  planeNormal(tolerance = 1e-7): Vector3 | undefined {
    const points = this.sample(17, 17).flat();
    const origin = points[0]!;
    let far = origin;
    for (const p of points) {
      if (p.distanceTo(origin) > far.distanceTo(origin)) {
        far = p;
      }
    }
    const axis = far.sub(origin);
    if (axis.length === 0) {
      return undefined;
    }
    let normal = Vector3.ZERO;
    for (const p of points) {
      const candidate = axis.cross(p.sub(origin));
      if (candidate.length > normal.length) {
        normal = candidate;
      }
    }
    if (normal.length <= tolerance * axis.length) {
      return undefined;
    }
    const unit = normal.normalized();
    if (points.some((p) => Math.abs(p.sub(origin).dot(unit)) > tolerance)) {
      return undefined;
    }
    const facing = this.normalAt(this.atU(0.37), this.atV(0.41));
    return facing.dot(unit) < 0 ? unit.negate() : unit;
  }

  isPlanar(tolerance = 1e-7): boolean {
    return this.planeNormal(tolerance) !== undefined;
  }

  /** Validates a trim rectangle, returning it clamped and in increasing order. */
  protected checkedRange(u0: number, u1: number, v0: number, v1: number): [number, number, number, number] {
    const { domainU: du, domainV: dv } = this;
    const a = clamp(Math.min(u0, u1), du.start, du.end);
    const b = clamp(Math.max(u0, u1), du.start, du.end);
    const c = clamp(Math.min(v0, v1), dv.start, dv.end);
    const d = clamp(Math.max(v0, v1), dv.start, dv.end);
    if (!(b > a) || !(d > c)) {
      throw new RangeError('Trim range must have a positive size inside the domain.');
    }
    return [a, b, c, d];
  }

  trim(u0: number, u1: number, v0: number, v1: number): ISurface {
    return new TrimmedSurface(this, ...this.checkedRange(u0, u1, v0, v1));
  }

  splitU(u: number): [ISurface, ISurface] {
    const { domainU: du, domainV: dv } = this;
    if (!(u > du.start && u < du.end)) {
      throw new RangeError('Split parameter must be strictly inside the domain.');
    }
    return [this.trim(du.start, u, dv.start, dv.end), this.trim(u, du.end, dv.start, dv.end)];
  }

  splitV(v: number): [ISurface, ISurface] {
    const { domainU: du, domainV: dv } = this;
    if (!(v > dv.start && v < dv.end)) {
      throw new RangeError('Split parameter must be strictly inside the domain.');
    }
    return [this.trim(du.start, du.end, dv.start, v), this.trim(du.start, du.end, v, dv.end)];
  }

  reversed(): ISurface {
    return new ReversedSurface(this);
  }

  transformed(transform: Transform3): ISurface {
    return new TransformedSurface(this, transform);
  }

  offset(distance: number): ISurface {
    return new OffsetSurface(this, distance);
  }

  toNurbs(tolerance = 1e-4): NurbsSurface {
    if (!nurbsApproximator) {
      throw new Error('NURBS support is not loaded.');
    }
    return nurbsApproximator(this, tolerance);
  }

  equals(other: ISurface, tolerance = 1e-7): boolean {
    const a = this.sample(9, 9);
    const b = other.sample(9, 9);
    return a.every((row, i) => row.every((p, j) => p.distanceTo(b[i]![j]!) <= tolerance));
  }
}

/** A curve on a surface at one fixed parameter. Its derivatives come from the surface's. */
export class IsoCurve extends CurveBase {
  readonly type = 'IsoCurve';
  readonly domain: CurveDomain;

  constructor(
    readonly surface: ISurface,
    /** The parameter that is fixed. */
    readonly fixed: 'u' | 'v',
    readonly value: number,
  ) {
    super();
    const { start, end } = fixed === 'u' ? surface.domainV : surface.domainU;
    this.domain = { start, end };
  }

  override get isPeriodic(): boolean {
    return this.fixed === 'u' ? this.surface.isPeriodicV : this.surface.isPeriodicU;
  }

  protected evaluate(t: number): Vector3 {
    return this.fixed === 'u' ? this.surface.pointAt(this.value, t) : this.surface.pointAt(t, this.value);
  }

  protected override evaluateDerivative(t: number, order: number): Vector3 {
    return this.fixed === 'u' ? this.surface.derivativeAt(this.value, t, 0, order) : this.surface.derivativeAt(t, this.value, order, 0);
  }

  override breakpoints(): number[] {
    return this.fixed === 'u' ? this.surface.breakpointsV() : this.surface.breakpointsU();
  }

  toJSON(): CurveJSON {
    return { type: this.type, surface: this.surface.toJSON(), fixed: this.fixed, value: this.value };
  }
}

/** A surface over a sub-rectangle of another's domain, with the same parameters. */
export class TrimmedSurface extends SurfaceBase {
  readonly type = 'TrimmedSurface';
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;

  constructor(
    readonly surface: ISurface,
    u0: number,
    u1: number,
    v0: number,
    v1: number,
  ) {
    super();
    this.domainU = { start: u0, end: u1 };
    this.domainV = { start: v0, end: v1 };
  }

  override get degreeU(): number | undefined {
    return this.surface.degreeU;
  }

  override get degreeV(): number | undefined {
    return this.surface.degreeV;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.surface.pointAt(u, v);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.surface.derivativeAt(u, v, orderU, orderV);
  }

  override breakpointsU(): number[] {
    return clipBreaks(this.surface.breakpointsU(), this.domainU.start, this.domainU.end);
  }

  override breakpointsV(): number[] {
    return clipBreaks(this.surface.breakpointsV(), this.domainV.start, this.domainV.end);
  }

  override trim(u0: number, u1: number, v0: number, v1: number): ISurface {
    return new TrimmedSurface(this.surface, ...this.checkedRange(u0, u1, v0, v1));
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      surface: this.surface.toJSON(),
      u: [this.domainU.start, this.domainU.end],
      v: [this.domainV.start, this.domainV.end],
    };
  }
}

/** The same geometry with `u` mirrored over its domain, so the normal is flipped. */
export class ReversedSurface extends SurfaceBase {
  readonly type = 'ReversedSurface';

  constructor(readonly surface: ISurface) {
    super();
  }

  get domainU(): SurfaceDomain {
    return this.surface.domainU;
  }

  get domainV(): SurfaceDomain {
    return this.surface.domainV;
  }

  override get degreeU(): number | undefined {
    return this.surface.degreeU;
  }

  override get degreeV(): number | undefined {
    return this.surface.degreeV;
  }

  override get isPeriodicU(): boolean {
    return this.surface.isPeriodicU;
  }

  override get isPeriodicV(): boolean {
    return this.surface.isPeriodicV;
  }

  private flip(u: number): number {
    return this.domainU.start + this.domainU.end - u;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.surface.pointAt(this.flip(u), v);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const d = this.surface.derivativeAt(this.flip(u), v, orderU, orderV);
    return orderU % 2 === 1 ? d.negate() : d;
  }

  override breakpointsU(): number[] {
    return this.surface.breakpointsU().map((b) => this.flip(b)).reverse();
  }

  override breakpointsV(): number[] {
    return this.surface.breakpointsV();
  }

  override reversed(): ISurface {
    return this.surface;
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, surface: this.surface.toJSON() };
  }
}

/** Another surface under an affine transform, with the same parameters. */
export class TransformedSurface extends SurfaceBase {
  readonly type = 'TransformedSurface';

  constructor(
    readonly surface: ISurface,
    readonly transform: Transform3,
  ) {
    super();
  }

  get domainU(): SurfaceDomain {
    return this.surface.domainU;
  }

  get domainV(): SurfaceDomain {
    return this.surface.domainV;
  }

  override get degreeU(): number | undefined {
    return this.surface.degreeU;
  }

  override get degreeV(): number | undefined {
    return this.surface.degreeV;
  }

  override get isPeriodicU(): boolean {
    return this.surface.isPeriodicU;
  }

  override get isPeriodicV(): boolean {
    return this.surface.isPeriodicV;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.transform.applyToPoint(this.surface.pointAt(u, v));
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.transform.applyToVector(this.surface.derivativeAt(u, v, orderU, orderV));
  }

  override breakpointsU(): number[] {
    return this.surface.breakpointsU();
  }

  override breakpointsV(): number[] {
    return this.surface.breakpointsV();
  }

  override transformed(transform: Transform3): ISurface {
    return new TransformedSurface(this.surface, transform.compose(this.transform));
  }

  /** Exact whenever the original's NURBS form is, since affine maps keep NURBS shapes. */
  override toNurbs(tolerance = 1e-4): NurbsSurface {
    return this.surface.toNurbs(tolerance).transformed(this.transform);
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      surface: this.surface.toJSON(),
      linear: [...this.transform.linear],
      translation: this.transform.translation.toArray(),
    };
  }
}

/**
 * Points at a fixed distance along another surface's normal, with the same parameters. It follows
 * the base exactly, so it folds over itself where the distance exceeds a radius of curvature.
 */
export class OffsetSurface extends SurfaceBase {
  readonly type = 'OffsetSurface';

  constructor(
    readonly surface: ISurface,
    readonly distance: number,
  ) {
    super();
  }

  get domainU(): SurfaceDomain {
    return this.surface.domainU;
  }

  get domainV(): SurfaceDomain {
    return this.surface.domainV;
  }

  override get isPeriodicU(): boolean {
    return this.surface.isPeriodicU;
  }

  override get isPeriodicV(): boolean {
    return this.surface.isPeriodicV;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.surface.pointAt(u, v).addScaled(this.surface.normalAt(u, v), this.distance);
  }

  override breakpointsU(): number[] {
    return this.surface.breakpointsU();
  }

  override breakpointsV(): number[] {
    return this.surface.breakpointsV();
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, surface: this.surface.toJSON(), distance: this.distance };
  }
}

/** Determinant of a transform's linear part. Negative for a reflection. */
export function linearDeterminant(t: Transform3): number {
  const m = t.linear;
  return (
    m[0]! * (m[4]! * m[8]! - m[5]! * m[7]!) -
    m[1]! * (m[3]! * m[8]! - m[5]! * m[6]!) +
    m[2]! * (m[3]! * m[7]! - m[4]! * m[6]!)
  );
}
