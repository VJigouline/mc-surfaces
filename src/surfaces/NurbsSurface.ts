import { NurbsCurve, Vector3, type ICurve, type Transform3 } from 'meercad-curves';
import { averagedKnots, basisDerivatives, basisRow, findSpan, uniformClampedKnots } from '../math/bspline.js';
import { binomial, solveVectorSystem } from '../math/numeric.js';
import type { ISurface, SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { registerNurbsApproximator, SurfaceBase } from './SurfaceBase.js';

function checkKnots(degree: number, count: number, knots: readonly number[], name: string): void {
  if (!Number.isInteger(degree) || degree < 1) {
    throw new RangeError('NURBS degree must be a positive integer.');
  }
  if (count < degree + 1) {
    throw new RangeError(`Degree ${degree} needs at least ${degree + 1} control points in ${name}.`);
  }
  if (knots.length !== count + degree + 1) {
    throw new RangeError(`Expected ${count + degree + 1} ${name} knots, got ${knots.length}.`);
  }
  for (let i = 0; i + 1 < knots.length; i++) {
    if (!(knots[i]! <= knots[i + 1]!)) {
      throw new RangeError('Knots must be non-decreasing.');
    }
  }
  if (!(knots[count]! > knots[degree]!)) {
    throw new RangeError('The knot vector leaves an empty domain.');
  }
}

function distinctKnots(knots: readonly number[], domain: SurfaceDomain): number[] {
  const result = [domain.start];
  for (const k of knots) {
    if (k > result[result.length - 1]! && k < domain.end) {
      result.push(k);
    }
  }
  result.push(domain.end);
  return result;
}

/**
 * Tensor-product non-uniform rational B-spline. Omit `weights` for a plain B-spline.
 * Control points and weights are indexed `[u][v]`. Algorithm numbers refer to The NURBS Book.
 */
export class NurbsSurface extends SurfaceBase {
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;
  readonly controlPoints: readonly (readonly Vector3[])[];
  readonly knotsU: readonly number[];
  readonly knotsV: readonly number[];
  readonly weights: readonly (readonly number[])[] | undefined;
  private readonly splineDegreeU: number;
  private readonly splineDegreeV: number;

  constructor(
    degreeU: number,
    degreeV: number,
    controlPoints: readonly (readonly Vector3[])[],
    knotsU: readonly number[],
    knotsV: readonly number[],
    weights?: readonly (readonly number[])[],
  ) {
    super();
    const countU = controlPoints.length;
    const countV = controlPoints[0]?.length ?? 0;
    if (controlPoints.some((row) => row.length !== countV)) {
      throw new RangeError('Every row of control points must have the same length.');
    }
    checkKnots(degreeU, countU, knotsU, 'u');
    checkKnots(degreeV, countV, knotsV, 'v');
    if (weights) {
      if (weights.length !== countU || weights.some((row) => row.length !== countV)) {
        throw new RangeError('There must be one weight per control point.');
      }
      if (weights.some((row) => row.some((w) => !(w > 0)))) {
        throw new RangeError('Weights must be positive.');
      }
    }
    this.splineDegreeU = degreeU;
    this.splineDegreeV = degreeV;
    this.controlPoints = controlPoints.map((row) => [...row]);
    this.knotsU = [...knotsU];
    this.knotsV = [...knotsV];
    const first = weights?.[0]?.[0];
    this.weights = weights && weights.some((row) => row.some((w) => w !== first)) ? weights.map((row) => [...row]) : undefined;
    this.domainU = { start: knotsU[degreeU]!, end: knotsU[countU]! };
    this.domainV = { start: knotsV[degreeV]!, end: knotsV[countV]! };
  }

  get type(): string {
    return this.isRational ? 'NurbsSurface' : 'BSplineSurface';
  }

  override get degreeU(): number {
    return this.splineDegreeU;
  }

  override get degreeV(): number {
    return this.splineDegreeV;
  }

  get isRational(): boolean {
    return this.weights !== undefined;
  }

  weightAt(i: number, j: number): number {
    return this.weights?.[i]?.[j] ?? 1;
  }

  /** All derivatives up to the given orders, indexed `[k][l]`, using eq. 4.20 for rational forms. */
  private derivativesAt(u: number, v: number, orderU: number, orderV: number): Vector3[][] {
    const p = this.splineDegreeU;
    const q = this.splineDegreeV;
    const countU = this.controlPoints.length;
    const countV = this.controlPoints[0]!.length;
    const spanU = findSpan(u, p, this.knotsU, countU);
    const spanV = findSpan(v, q, this.knotsV, countV);
    const basisU = basisDerivatives(u, spanU, p, this.knotsU, orderU);
    const basisV = basisDerivatives(v, spanV, q, this.knotsV, orderV);
    const numerator: Vector3[][] = [];
    const weight: number[][] = [];
    for (let k = 0; k <= orderU; k++) {
      numerator.push([]);
      weight.push([]);
      for (let l = 0; l <= orderV; l++) {
        let [x, y, z, w] = [0, 0, 0, 0];
        for (let a = 0; a <= p; a++) {
          const row = this.controlPoints[spanU - p + a]!;
          const nu = basisU[k]![a]!;
          if (nu === 0) {
            continue;
          }
          for (let b = 0; b <= q; b++) {
            const pt = row[spanV - q + b]!;
            const n = nu * basisV[l]![b]!;
            const wt = this.weights ? this.weights[spanU - p + a]![spanV - q + b]! : 1;
            x += n * wt * pt.x;
            y += n * wt * pt.y;
            z += n * wt * pt.z;
            w += n * wt;
          }
        }
        numerator[k]!.push(new Vector3(x, y, z));
        weight[k]!.push(w);
      }
    }
    if (!this.weights) {
      return numerator;
    }
    const result: Vector3[][] = [];
    for (let k = 0; k <= orderU; k++) {
      result.push([]);
      for (let l = 0; l <= orderV; l++) {
        let s = numerator[k]![l]!;
        for (let i = 1; i <= k; i++) {
          s = s.addScaled(result[k - i]![l]!, -binomial(k, i) * weight[i]![0]!);
        }
        for (let j = 1; j <= l; j++) {
          s = s.addScaled(result[k]![l - j]!, -binomial(l, j) * weight[0]![j]!);
        }
        for (let i = 1; i <= k; i++) {
          for (let j = 1; j <= l; j++) {
            s = s.addScaled(result[k - i]![l - j]!, -binomial(k, i) * binomial(l, j) * weight[i]![j]!);
          }
        }
        result[k]!.push(s.scale(1 / weight[0]![0]!));
      }
    }
    return result;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.derivativesAt(u, v, 0, 0)[0]![0]!;
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.derivativesAt(u, v, orderU, orderV)[orderU]![orderV]!;
  }

  override breakpointsU(): number[] {
    return distinctKnots(this.knotsU, this.domainU);
  }

  override breakpointsV(): number[] {
    return distinctKnots(this.knotsV, this.domainV);
  }

  /** The exact curve at fixed `u`, by blending rows in homogeneous form. */
  override isoCurveU(u: number): ICurve {
    const t = this.normalizeU(u);
    const p = this.splineDegreeU;
    const span = findSpan(t, p, this.knotsU, this.controlPoints.length);
    const basis = basisDerivatives(t, span, p, this.knotsU, 0)[0]!;
    const points: Vector3[] = [];
    const weights: number[] = [];
    for (let j = 0; j < this.controlPoints[0]!.length; j++) {
      let sum = Vector3.ZERO;
      let w = 0;
      for (let a = 0; a <= p; a++) {
        const i = span - p + a;
        const wt = basis[a]! * this.weightAt(i, j);
        sum = sum.addScaled(this.controlPoints[i]![j]!, wt);
        w += wt;
      }
      points.push(sum.scale(1 / w));
      weights.push(w);
    }
    return new NurbsCurve(this.splineDegreeV, points, this.knotsV, this.weights ? weights : undefined);
  }

  /** The exact curve at fixed `v`, by blending columns in homogeneous form. */
  override isoCurveV(v: number): ICurve {
    const t = this.normalizeV(v);
    const q = this.splineDegreeV;
    const span = findSpan(t, q, this.knotsV, this.controlPoints[0]!.length);
    const basis = basisDerivatives(t, span, q, this.knotsV, 0)[0]!;
    const points: Vector3[] = [];
    const weights: number[] = [];
    for (let i = 0; i < this.controlPoints.length; i++) {
      let sum = Vector3.ZERO;
      let w = 0;
      for (let b = 0; b <= q; b++) {
        const j = span - q + b;
        const wt = basis[b]! * this.weightAt(i, j);
        sum = sum.addScaled(this.controlPoints[i]![j]!, wt);
        w += wt;
      }
      points.push(sum.scale(1 / w));
      weights.push(w);
    }
    return new NurbsCurve(this.splineDegreeU, points, this.knotsU, this.weights ? weights : undefined);
  }

  /** Applies a curve operation to every line of control points running in `u`. */
  private mapAlongU(operation: (curve: NurbsCurve) => NurbsCurve): NurbsSurface {
    const countV = this.controlPoints[0]!.length;
    const columns: NurbsCurve[] = [];
    const columnWeights: number[][] = [];
    for (let j = 0; j < countV; j++) {
      const w = this.controlPoints.map((_, i) => this.weightAt(i, j));
      const result = operation(new NurbsCurve(this.splineDegreeU, this.controlPoints.map((row) => row[j]!), this.knotsU, w));
      columns.push(result);
      // A curve drops weights that are all equal, so restore that common value.
      columnWeights.push(result.controlPoints.map((_, i) => result.weights?.[i] ?? w[0]!));
    }
    const points = columns[0]!.controlPoints.map((_, i) => columns.map((c) => c.controlPoints[i]!));
    const weights = columns[0]!.controlPoints.map((_, i) => columnWeights.map((c) => c[i]!));
    return new NurbsSurface(this.splineDegreeU, this.splineDegreeV, points, columns[0]!.knots, this.knotsV, this.weights ? weights : undefined);
  }

  /** Applies a curve operation to every line of control points running in `v`. */
  private mapAlongV(operation: (curve: NurbsCurve) => NurbsCurve): NurbsSurface {
    const rows = this.controlPoints.map((row, i) => {
      const w = row.map((_, j) => this.weightAt(i, j));
      const result = operation(new NurbsCurve(this.splineDegreeV, row, this.knotsV, w));
      return { result, weights: result.controlPoints.map((_, j) => result.weights?.[j] ?? w[0]!) };
    });
    return new NurbsSurface(
      this.splineDegreeU,
      this.splineDegreeV,
      rows.map((r) => r.result.controlPoints),
      this.knotsU,
      rows[0]!.result.knots,
      this.weights ? rows.map((r) => r.weights) : undefined,
    );
  }

  /** The same surface with the knot `u` inserted `times` more. */
  insertKnotU(u: number, times = 1): NurbsSurface {
    return this.mapAlongU((curve) => curve.insertKnot(u, times));
  }

  /** The same surface with the knot `v` inserted `times` more. */
  insertKnotV(v: number, times = 1): NurbsSurface {
    return this.mapAlongV((curve) => curve.insertKnot(v, times));
  }

  override trim(u0: number, u1: number, v0: number, v1: number): ISurface {
    const [a, b, c, d] = this.checkedRange(u0, u1, v0, v1);
    return this.mapAlongU((curve) => curve.trim(a, b)).mapAlongV((curve) => curve.trim(c, d));
  }

  override reversed(): ISurface {
    const { start, end } = this.domainU;
    const knots = [...this.knotsU].reverse().map((k) => start + end - k);
    const points = [...this.controlPoints].reverse();
    const weights = this.weights ? [...this.weights].reverse() : undefined;
    return new NurbsSurface(this.splineDegreeU, this.splineDegreeV, points, knots, this.knotsV, weights);
  }

  /** Any affine map keeps a NURBS surface, by moving its control points. */
  override transformed(transform: Transform3): NurbsSurface {
    return this.mapControlPoints((p) => transform.applyToPoint(p));
  }

  /** The same knots and weights, with each control point replaced. The map must be affine to keep the shape's meaning. */
  mapControlPoints(map: (p: Vector3) => Vector3): NurbsSurface {
    return new NurbsSurface(
      this.splineDegreeU,
      this.splineDegreeV,
      this.controlPoints.map((row) => row.map(map)),
      this.knotsU,
      this.knotsV,
      this.weights,
    );
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    return this;
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      degreeU: this.splineDegreeU,
      degreeV: this.splineDegreeV,
      controlPoints: this.controlPoints.map((row) => row.map((p) => p.toArray())),
      knotsU: [...this.knotsU],
      knotsV: [...this.knotsV],
      ...(this.weights ? { weights: this.weights.map((row) => [...row]) } : {}),
    };
  }
}

/** A non-rational `NurbsSurface`, with two extra constructors. */
export class BSplineSurface extends NurbsSurface {
  constructor(
    degreeU: number,
    degreeV: number,
    controlPoints: readonly (readonly Vector3[])[],
    knotsU: readonly number[],
    knotsV: readonly number[],
  ) {
    super(degreeU, degreeV, controlPoints, knotsU, knotsV);
  }

  /** Clamped, with even knots over [0, 1] in both directions. */
  static uniform(degreeU: number, degreeV: number, controlPoints: readonly (readonly Vector3[])[]): BSplineSurface {
    const countU = controlPoints.length;
    const countV = controlPoints[0]?.length ?? 0;
    if (countU < degreeU + 1 || countV < degreeV + 1) {
      throw new RangeError(`Degrees (${degreeU}, ${degreeV}) need at least ${degreeU + 1} by ${degreeV + 1} control points.`);
    }
    return new BSplineSurface(degreeU, degreeV, controlPoints, uniformClampedKnots(degreeU, countU), uniformClampedKnots(degreeV, countV));
  }

  /**
   * Global interpolation through a grid of points, over [0, 1] in both directions, with averaged
   * chord-length parameters. Algorithm A9.4. Cost grows up to cubically with each grid side.
   */
  static interpolate(points: readonly (readonly Vector3[])[], degreeU = 3, degreeV = 3): BSplineSurface {
    const countU = points.length;
    const countV = points[0]?.length ?? 0;
    if (countU < degreeU + 1 || countV < degreeV + 1) {
      throw new RangeError(`Degrees (${degreeU}, ${degreeV}) need at least ${degreeU + 1} by ${degreeV + 1} points.`);
    }
    if (points.some((row) => row.length !== countV)) {
      throw new RangeError('Every row of points must have the same length.');
    }
    const line = (i: number, across: boolean) => (across ? points.map((row) => row[i]!) : points[i]!);
    const averaged = (lines: number, length: number, across: boolean): number[] => {
      const sums = new Array<number>(length).fill(0);
      let used = 0;
      for (let k = 0; k < lines; k++) {
        const pts = line(k, across);
        const chords = pts.slice(1).map((p, i) => p.distanceTo(pts[i]!));
        const total = chords.reduce((a, b) => a + b, 0);
        if (total === 0) {
          continue;
        }
        used++;
        let run = 0;
        chords.forEach((c, i) => {
          run += c;
          sums[i + 1]! += run / total;
        });
      }
      if (used === 0) {
        throw new RangeError('Interpolation points must not all coincide.');
      }
      const params = sums.map((s) => s / used);
      params[length - 1] = 1;
      for (let i = 0; i + 1 < length; i++) {
        if (!(params[i + 1]! > params[i]!)) {
          throw new RangeError('Interpolation points must not repeat consecutively.');
        }
      }
      return params;
    };
    const paramsU = averaged(countV, countU, true);
    const paramsV = averaged(countU, countV, false);
    const knotsU = averagedKnots(paramsU, degreeU);
    const knotsV = averagedKnots(paramsV, degreeV);
    const matrixU = paramsU.map((u) => basisRow(u, degreeU, knotsU, countU));
    const matrixV = paramsV.map((v) => basisRow(v, degreeV, knotsV, countV));
    const stage = Array.from({ length: countV }, (_, j) => solveVectorSystem(matrixU, points.map((row) => row[j]!)));
    const control = Array.from({ length: countU }, (_, i) => solveVectorSystem(matrixV, stage.map((col) => col[i]!)));
    return new BSplineSurface(degreeU, degreeV, control, knotsU, knotsV);
  }
}

/**
 * A tensor-product Bezier patch over [0, 1] in both directions, of any degree.
 * Pass `weights` for the rational form, reported as `RationalBezierSurface`.
 */
export class BezierSurface extends NurbsSurface {
  constructor(controlPoints: readonly (readonly Vector3[])[], weights?: readonly (readonly number[])[]) {
    const countU = controlPoints.length;
    const countV = controlPoints[0]?.length ?? 0;
    if (countU < 2 || countV < 2) {
      throw new RangeError('A Bezier surface needs at least 2 by 2 control points.');
    }
    const clamped = (count: number) => [...new Array<number>(count).fill(0), ...new Array<number>(count).fill(1)];
    super(countU - 1, countV - 1, controlPoints, clamped(countU), clamped(countV), weights);
  }

  override get type(): string {
    return this.isRational ? 'RationalBezierSurface' : 'BezierSurface';
  }

  override transformed(transform: Transform3): BezierSurface {
    return new BezierSurface(this.controlPoints.map((row) => row.map((p) => transform.applyToPoint(p))), this.weights);
  }

  override toJSON(): SurfaceJSON {
    return {
      type: this.type,
      controlPoints: this.controlPoints.map((row) => row.map((p) => p.toArray())),
      ...(this.weights ? { weights: this.weights.map((row) => [...row]) } : {}),
    };
  }
}

/** Control points and weights of a degree 2 arc, as `(cos, sin, weight)` of unit-circle points. */
function arcTemplate(start: number, sweep: number): { rows: [number, number, number][]; knots: number[] } {
  const pieces = Math.max(1, Math.ceil(sweep / (Math.PI / 2) - 1e-12));
  const step = sweep / pieces;
  const w = Math.cos(step / 2);
  const rows: [number, number, number][] = [];
  const knots = [start, start, start];
  for (let k = 0; k < pieces; k++) {
    const a = start + k * step;
    rows.push([Math.cos(a), Math.sin(a), 1]);
    const mid = a + step / 2;
    rows.push([Math.cos(mid) / w, Math.sin(mid) / w, w]);
    if (k > 0) {
      knots.push(a, a);
    }
  }
  const end = start + sweep;
  rows.push([Math.cos(end), Math.sin(end), 1]);
  knots.push(end, end, end);
  return { rows, knots };
}

/**
 * Exact NURBS surface of revolution, `u` being the angle and `v` the profile's parameter.
 * Algorithm A8.1, generalized to any start angle. `start` is the rotation of the profile as given.
 */
export function revolveNurbs(profile: NurbsCurve, axisOrigin: Vector3, axisDirection: Vector3, start: number, end: number): NurbsSurface {
  const axis = axisDirection.normalized();
  const { rows, knots } = arcTemplate(start, end - start);
  const points: Vector3[][] = [];
  const weights: number[][] = [];
  for (const [c, s, w] of rows) {
    const row: Vector3[] = [];
    const rowWeights: number[] = [];
    profile.controlPoints.forEach((p, j) => {
      const relative = p.sub(axisOrigin);
      const onAxis = axisOrigin.addScaled(axis, relative.dot(axis));
      const radial = p.sub(onAxis);
      row.push(onAxis.addScaled(radial, c).addScaled(axis.cross(radial), s));
      rowWeights.push(w * (profile.weights?.[j] ?? 1));
    });
    points.push(row);
    weights.push(rowWeights);
  }
  return new NurbsSurface(2, profile.degree, points, knots, profile.knots, weights);
}

const MAX_APPROXIMATION_CELLS = 256;

/**
 * Bicubic Hermite patches from the surface's point, tangents and twist at a grid of parameters,
 * written as Bezier spans. C1 inside smooth spans; the grid is refined until the check samples pass.
 */
function approximate(surface: ISurface, tolerance: number): NurbsSurface {
  if (!(tolerance > 0)) {
    throw new RangeError('NURBS tolerance must be positive.');
  }
  const refine = (breaks: readonly number[], pieces: number) => {
    const result: number[] = [];
    for (let i = 0; i + 1 < breaks.length; i++) {
      for (let k = 0; k < pieces; k++) {
        result.push(breaks[i]! + ((breaks[i + 1]! - breaks[i]!) * k) / pieces);
      }
    }
    result.push(breaks[breaks.length - 1]!);
    return result;
  };
  const breaksU = surface.breakpointsU();
  const breaksV = surface.breakpointsV();
  let best: NurbsSurface | undefined;
  for (let pieces = 2; ; pieces *= 2) {
    const us = refine(breaksU, pieces);
    const vs = refine(breaksV, pieces);
    const candidate = bicubicNet(surface, us, vs, breaksU, breaksV);
    best = candidate;
    const tooLarge = us.length * 2 > MAX_APPROXIMATION_CELLS || vs.length * 2 > MAX_APPROXIMATION_CELLS;
    if (tooLarge || withinTolerance(surface, candidate, us, vs, tolerance)) {
      return best;
    }
  }
}

function withinTolerance(surface: ISurface, nurbs: NurbsSurface, us: readonly number[], vs: readonly number[], tolerance: number): boolean {
  const fractions = [0.25, 0.5, 0.75];
  for (let i = 0; i + 1 < us.length; i++) {
    for (let j = 0; j + 1 < vs.length; j++) {
      for (const a of fractions) {
        for (const b of fractions) {
          const u = us[i]! + (us[i + 1]! - us[i]!) * a;
          const v = vs[j]! + (vs[j + 1]! - vs[j]!) * b;
          if (surface.pointAt(u, v).distanceTo(nurbs.pointAt(u, v)) > tolerance) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

function bicubicNet(
  surface: ISurface,
  us: readonly number[],
  vs: readonly number[],
  breaksU: readonly number[],
  breaksV: readonly number[],
): NurbsSurface {
  const spanU = surface.domainU.end - surface.domainU.start;
  const spanV = surface.domainV.end - surface.domainV.start;
  // One-sided data: a derivative taken just before an interior breakpoint belongs to the span on its left.
  const side = (t: number, s: number, breaks: readonly number[], span: number) =>
    s < 0 && breaks.slice(1, -1).some((b) => Math.abs(b - t) < 1e-12 * span) ? t - 1e-11 * span : t;
  const countU = 3 * (us.length - 1) + 1;
  const countV = 3 * (vs.length - 1) + 1;
  const net: Vector3[][] = Array.from({ length: countU }, () => new Array<Vector3>(countV).fill(Vector3.ZERO));
  us.forEach((u, i) => {
    vs.forEach((v, j) => {
      const p = surface.pointAt(u, v);
      for (const su of [-1, 1]) {
        const a = 3 * i + su;
        if (a < 0 || a >= countU) {
          continue;
        }
        const du = su < 0 ? (u - us[i - 1]!) / 3 : (us[i + 1]! - u) / 3;
        const atU = side(u, su, breaksU, spanU);
        const tangentU = surface.derivativeAt(atU, v, 1, 0);
        net[a]![3 * j] = p.addScaled(tangentU, su * du);
        for (const sv of [-1, 1]) {
          const b = 3 * j + sv;
          if (b < 0 || b >= countV) {
            continue;
          }
          const dv = sv < 0 ? (v - vs[j - 1]!) / 3 : (vs[j + 1]! - v) / 3;
          const atV = side(v, sv, breaksV, spanV);
          const tangentV = surface.derivativeAt(atU, atV, 0, 1);
          const twist = surface.derivativeAt(atU, atV, 1, 1);
          const tangentAtSide = surface.derivativeAt(atU, atV, 1, 0);
          net[a]![b] = p.addScaled(tangentAtSide, su * du).addScaled(tangentV, sv * dv).addScaled(twist, su * du * sv * dv);
        }
      }
      for (const sv of [-1, 1]) {
        const b = 3 * j + sv;
        if (b < 0 || b >= countV) {
          continue;
        }
        const dv = sv < 0 ? (v - vs[j - 1]!) / 3 : (vs[j + 1]! - v) / 3;
        net[3 * i]![b] = p.addScaled(surface.derivativeAt(u, side(v, sv, breaksV, spanV), 0, 1), sv * dv);
      }
      net[3 * i]![3 * j] = p;
    });
  });
  const knots = (values: readonly number[]) => {
    const result = [values[0]!, values[0]!, values[0]!, values[0]!];
    for (let i = 1; i + 1 < values.length; i++) {
      result.push(values[i]!, values[i]!, values[i]!);
    }
    const last = values[values.length - 1]!;
    result.push(last, last, last, last);
    return result;
  };
  return new NurbsSurface(3, 3, net, knots(us), knots(vs));
}

registerNurbsApproximator(approximate);
