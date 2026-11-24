import type { ICurve, Transform3, Vector3 } from 'mc-curves';
import { Vector3 as V3 } from 'mc-curves';
import { CardinalSplines } from '../math/cardinal.js';
import { dCos, dSin, uniqueSorted } from '../math/numeric.js';
import { isFullTurn } from './FrameSurface.js';
import type { ISurface, SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { NurbsSurface, revolveNurbs } from './NurbsSurface.js';
import { linearDeterminant, SurfaceBase, TransformedSurface } from './SurfaceBase.js';

/** A curve read over [0, 1], whatever its own domain. */
class UnitCurve {
  private readonly start: number;
  private readonly span: number;

  constructor(readonly curve: ICurve) {
    this.start = curve.domain.start;
    this.span = curve.domain.end - curve.domain.start;
  }

  at(t: number, order: number): Vector3 {
    const c = this.start + this.span * t;
    return order === 0 ? this.curve.pointAt(c) : this.curve.derivativeAt(c, order).scale(this.span ** order);
  }

  breaks(): number[] {
    return this.curve.breakpoints().map((b) => (b - this.start) / this.span);
  }
}

function mergedBreaks(curves: readonly UnitCurve[], extra: readonly number[] = []): number[] {
  return uniqueSorted([0, 1, ...extra, ...curves.flatMap((c) => c.breaks())], 1e-10);
}

/** `C(u) + v D`: a curve moved along a direction. */
export class ExtrusionSurface extends SurfaceBase {
  readonly type = 'ExtrusionSurface';
  readonly domainV: SurfaceDomain;

  constructor(
    readonly curve: ICurve,
    readonly direction: Vector3,
    vStart = 0,
    vEnd = 1,
  ) {
    super();
    if (direction.lengthSquared === 0) {
      throw new RangeError('Extrusion direction must not be zero.');
    }
    if (!(vEnd > vStart)) {
      throw new RangeError('Surface domain must be of positive size.');
    }
    this.domainV = { start: vStart, end: vEnd };
  }

  get domainU(): SurfaceDomain {
    return this.curve.domain;
  }

  override get degreeU(): number | undefined {
    return this.curve.degree;
  }

  override get degreeV(): number {
    return 1;
  }

  override get isPeriodicU(): boolean {
    return this.curve.isPeriodic;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.curve.pointAt(u).addScaled(this.direction, v);
  }

  protected override evaluateDerivative(u: number, _v: number, orderU: number, orderV: number): Vector3 {
    if (orderV === 0) {
      return this.curve.derivativeAt(u, orderU);
    }
    return orderU === 0 && orderV === 1 ? this.direction : V3.ZERO;
  }

  override breakpointsU(): number[] {
    return this.curve.breakpoints();
  }

  override transformed(transform: Transform3): ExtrusionSurface {
    return new ExtrusionSurface(this.curve.transformed(transform), transform.applyToVector(this.direction), this.domainV.start, this.domainV.end);
  }

  /** Exact where the curve's NURBS form is, though the curve's parameter may change. */
  override toNurbs(tolerance = 1e-4): NurbsSurface {
    const c = this.curve.toNurbs(tolerance);
    const { start, end } = this.domainV;
    const points = c.controlPoints.map((p) => [p.addScaled(this.direction, start), p.addScaled(this.direction, end)]);
    const weights = c.controlPoints.map((_, i) => {
      const w = c.weights?.[i] ?? 1;
      return [w, w];
    });
    return new NurbsSurface(c.degree, 1, points, c.knots, [start, start, end, end], weights);
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, curve: this.curve.toJSON(), direction: this.direction.toArray(), v: [this.domainV.start, this.domainV.end] };
  }
}

/**
 * A curve turned about an axis, right-handed. `u` is the angle in radians, with `u = 0` at the
 * profile as given, and `v` is the curve's parameter.
 */
export class RevolutionSurface extends SurfaceBase {
  readonly type = 'RevolutionSurface';
  readonly domainU: SurfaceDomain;
  readonly axisDirection: Vector3;

  constructor(
    readonly curve: ICurve,
    readonly axisOrigin: Vector3,
    axisDirection: Vector3,
    uStart = 0,
    uEnd = 2 * Math.PI,
  ) {
    super();
    if (axisDirection.lengthSquared === 0) {
      throw new RangeError('Axis direction must not be zero.');
    }
    if (!(uEnd > uStart)) {
      throw new RangeError('Surface domain must be of positive size.');
    }
    this.axisDirection = axisDirection.normalized();
    this.domainU = { start: uStart, end: uEnd };
  }

  get domainV(): SurfaceDomain {
    return this.curve.domain;
  }

  override get degreeV(): number | undefined {
    return this.curve.degree;
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override get isPeriodicV(): boolean {
    return this.curve.isPeriodic;
  }

  private rotate(radial: Vector3, u: number, order: number): Vector3 {
    return radial.scale(dCos(u, order)).addScaled(this.axisDirection.cross(radial), dSin(u, order));
  }

  protected evaluate(u: number, v: number): Vector3 {
    const relative = this.curve.pointAt(v).sub(this.axisOrigin);
    const axial = this.axisDirection.scale(relative.dot(this.axisDirection));
    return this.axisOrigin.add(axial).add(this.rotate(relative.sub(axial), u, 0));
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const d = orderV === 0 ? this.curve.pointAt(v).sub(this.axisOrigin) : this.curve.derivativeAt(v, orderV);
    const axial = this.axisDirection.scale(d.dot(this.axisDirection));
    const radial = d.sub(axial);
    return orderU === 0 ? axial.add(this.rotate(radial, u, 0)) : this.rotate(radial, u, orderU);
  }

  override breakpointsV(): number[] {
    return this.curve.breakpoints();
  }

  /** A rotation stays a rotation only under a transform that keeps angles and handedness. */
  override transformed(transform: Transform3): ISurface {
    if (transform.similarityScale() === undefined || linearDeterminant(transform) < 0) {
      return new TransformedSurface(this, transform);
    }
    return new RevolutionSurface(
      this.curve.transformed(transform),
      transform.applyToPoint(this.axisOrigin),
      transform.applyToVector(this.axisDirection),
      this.domainU.start,
      this.domainU.end,
    );
  }

  /** Exact where the curve's NURBS form is. */
  override toNurbs(tolerance = 1e-4): NurbsSurface {
    return revolveNurbs(this.curve.toNurbs(tolerance), this.axisOrigin, this.axisDirection, this.domainU.start, this.domainU.end);
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      curve: this.curve.toJSON(),
      axisOrigin: this.axisOrigin.toArray(),
      axisDirection: this.axisDirection.toArray(),
      u: [this.domainU.start, this.domainU.end],
    };
  }
}

const UNIT: SurfaceDomain = { start: 0, end: 1 };

/** Straight lines between two curves: `(1 - v) A(u) + v B(u)`, with both read over `u` in [0, 1]. */
export class RuledSurface extends SurfaceBase {
  readonly type = 'RuledSurface';
  readonly domainU = UNIT;
  readonly domainV = UNIT;
  private readonly a: UnitCurve;
  private readonly b: UnitCurve;

  constructor(
    readonly curveA: ICurve,
    readonly curveB: ICurve,
  ) {
    super();
    this.a = new UnitCurve(curveA);
    this.b = new UnitCurve(curveB);
  }

  override get degreeV(): number {
    return 1;
  }

  override get isPeriodicU(): boolean {
    return this.curveA.isPeriodic && this.curveB.isPeriodic;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.a.at(u, 0).lerp(this.b.at(u, 0), v);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    if (orderV > 1) {
      return V3.ZERO;
    }
    const a = this.a.at(u, orderU);
    const b = this.b.at(u, orderU);
    return orderV === 1 ? b.sub(a) : a.lerp(b, v);
  }

  override breakpointsU(): number[] {
    return mergedBreaks([this.a, this.b]);
  }

  override transformed(transform: Transform3): RuledSurface {
    return new RuledSurface(this.curveA.transformed(transform), this.curveB.transformed(transform));
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, curveA: this.curveA.toJSON(), curveB: this.curveB.toJSON() };
  }
}

export interface LoftOptions {
  /** Where each curve lies in `v`, increasing. For a closed loft, one more value for the return to the first curve. Default: by the average spacing of the curves, over [0, 1]. */
  readonly parameters?: readonly number[];
  /** Continue smoothly from the last curve back to the first. */
  readonly closed?: boolean;
}

/**
 * A C2 surface through a series of curves, all read over `u` in [0, 1] and running the same way.
 * Across them, `v` follows cardinal cubic splines. Two curves give a ruled surface.
 */
export class LoftSurface extends SurfaceBase {
  readonly type = 'LoftSurface';
  readonly domainU = UNIT;
  readonly domainV: SurfaceDomain;
  readonly closed: boolean;
  readonly parameters: readonly number[];
  private readonly sections: readonly UnitCurve[];
  private readonly blend: CardinalSplines;

  constructor(
    readonly curves: readonly ICurve[],
    options: LoftOptions = {},
  ) {
    super();
    if (curves.length < 2) {
      throw new RangeError('A loft needs at least two curves.');
    }
    this.closed = options.closed ?? false;
    this.sections = curves.map((c) => new UnitCurve(c));
    this.parameters = options.parameters ?? this.spacedParameters();
    if (this.parameters.length !== curves.length + (this.closed ? 1 : 0)) {
      throw new RangeError('There must be one parameter per curve, and one more if the loft is closed.');
    }
    this.blend = new CardinalSplines(this.parameters, this.closed);
    this.domainV = { start: this.parameters[0]!, end: this.parameters[this.parameters.length - 1]! };
  }

  /** Spacing in proportion to the average distance between neighbouring curves. */
  private spacedParameters(): number[] {
    const gaps: number[] = [];
    const count = this.sections.length + (this.closed ? 0 : -1);
    for (let i = 0; i < count; i++) {
      const [a, b] = [this.sections[i]!, this.sections[(i + 1) % this.sections.length]!];
      let sum = 0;
      for (let k = 0; k <= 8; k++) {
        sum += a.at(k / 8, 0).distanceTo(b.at(k / 8, 0));
      }
      gaps.push(sum / 9);
    }
    const total = gaps.reduce((s, g) => s + g, 0);
    if (gaps.some((g) => !(g > 0))) {
      throw new RangeError('Neighbouring loft curves must not coincide.');
    }
    const parameters = [0];
    gaps.forEach((g) => parameters.push(parameters[parameters.length - 1]! + g / total));
    parameters[parameters.length - 1] = 1;
    return parameters;
  }

  override get isPeriodicU(): boolean {
    return this.curves.every((c) => c.isPeriodic);
  }

  override get isPeriodicV(): boolean {
    return this.closed;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.evaluateDerivative(u, v, 0, 0);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const weights = this.blend.weights(v, orderV);
    let sum = V3.ZERO;
    this.sections.forEach((section, i) => {
      if (weights[i] !== 0) {
        sum = sum.addScaled(section.at(u, orderU), weights[i]!);
      }
    });
    return sum;
  }

  override breakpointsU(): number[] {
    return mergedBreaks(this.sections);
  }

  override breakpointsV(): number[] {
    return [...this.parameters];
  }

  override transformed(transform: Transform3): LoftSurface {
    return new LoftSurface(this.curves.map((c) => c.transformed(transform)), { parameters: this.parameters, closed: this.closed });
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, curves: this.curves.map((c) => c.toJSON()), parameters: [...this.parameters], closed: this.closed };
  }
}

function linearBlend(t: number, order: number, atStart: boolean): number {
  if (order === 0) {
    return atStart ? 1 - t : t;
  }
  return order === 1 ? (atStart ? -1 : 1) : 0;
}

/**
 * Bilinearly blended Coons patch filling four boundary curves, over [0, 1] in both directions.
 * `bottom` and `top` run in `u`; `left` and `right` run in `v`. They must meet at the corners.
 */
export class CoonsSurface extends SurfaceBase {
  readonly type = 'CoonsSurface';
  readonly domainU = UNIT;
  readonly domainV = UNIT;
  private readonly edges: readonly [UnitCurve, UnitCurve, UnitCurve, UnitCurve];
  private readonly corners: readonly [Vector3, Vector3, Vector3, Vector3];

  constructor(
    readonly bottom: ICurve,
    readonly top: ICurve,
    readonly left: ICurve,
    readonly right: ICurve,
    tolerance = 1e-7,
  ) {
    super();
    if (
      bottom.startPoint.distanceTo(left.startPoint) > tolerance ||
      bottom.endPoint.distanceTo(right.startPoint) > tolerance ||
      top.startPoint.distanceTo(left.endPoint) > tolerance ||
      top.endPoint.distanceTo(right.endPoint) > tolerance
    ) {
      throw new RangeError('Coons boundary curves must meet at the corners.');
    }
    this.edges = [new UnitCurve(bottom), new UnitCurve(top), new UnitCurve(left), new UnitCurve(right)];
    this.corners = [bottom.startPoint, bottom.endPoint, top.startPoint, top.endPoint];
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.evaluateDerivative(u, v, 0, 0);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const [bottom, top, left, right] = this.edges;
    const [p00, p10, p01, p11] = this.corners;
    const fu0 = linearBlend(u, orderU, true);
    const fu1 = linearBlend(u, orderU, false);
    const fv0 = linearBlend(v, orderV, true);
    const fv1 = linearBlend(v, orderV, false);
    return bottom.at(u, orderU).scale(fv0)
      .addScaled(top.at(u, orderU), fv1)
      .addScaled(left.at(v, orderV), fu0)
      .addScaled(right.at(v, orderV), fu1)
      .addScaled(p00, -fu0 * fv0)
      .addScaled(p10, -fu1 * fv0)
      .addScaled(p01, -fu0 * fv1)
      .addScaled(p11, -fu1 * fv1);
  }

  override breakpointsU(): number[] {
    return mergedBreaks([this.edges[0], this.edges[1]]);
  }

  override breakpointsV(): number[] {
    return mergedBreaks([this.edges[2], this.edges[3]]);
  }

  override transformed(transform: Transform3): CoonsSurface {
    const t = (c: ICurve) => c.transformed(transform);
    return new CoonsSurface(t(this.bottom), t(this.top), t(this.left), t(this.right));
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, bottom: this.bottom.toJSON(), top: this.top.toJSON(), left: this.left.toJSON(), right: this.right.toJSON() };
  }
}

export interface GordonOptions {
  /** Where each `v` curve crosses the `u` curves, from 0 to 1. Default: found from the first `u` curve. */
  readonly parametersU?: readonly number[];
  /** Where each `u` curve crosses the `v` curves, from 0 to 1. Default: found from the first `v` curve. */
  readonly parametersV?: readonly number[];
  /** Largest allowed gap where two curves cross. Default 1e-6 of the network's size. */
  readonly tolerance?: number;
}

/**
 * A surface through a network of crossing curves, over [0, 1] in both directions: the sum of a loft
 * through each family minus the surface through their crossing points.
 * Each `u` curve must cross each `v` curve at the same normalized parameters as its neighbours.
 */
export class GordonSurface extends SurfaceBase {
  readonly type = 'GordonSurface';
  readonly domainU = UNIT;
  readonly domainV = UNIT;
  readonly parametersU: readonly number[];
  readonly parametersV: readonly number[];
  private readonly families: { readonly u: readonly UnitCurve[]; readonly v: readonly UnitCurve[] };
  private readonly blendU: CardinalSplines;
  private readonly blendV: CardinalSplines;
  private readonly crossings: readonly (readonly Vector3[])[];

  constructor(
    readonly curvesU: readonly ICurve[],
    readonly curvesV: readonly ICurve[],
    options: GordonOptions = {},
  ) {
    super();
    if (curvesU.length < 2 || curvesV.length < 2) {
      throw new RangeError('A Gordon surface needs at least two curves in each family.');
    }
    const u = curvesU.map((c) => new UnitCurve(c));
    const v = curvesV.map((c) => new UnitCurve(c));
    const located = (on: ICurve, others: readonly ICurve[]) =>
      others.map((c, i) => {
        if (i === 0) {
          return 0;
        }
        if (i === others.length - 1) {
          return 1;
        }
        return (on.closestPoint(c.startPoint).parameter - on.domain.start) / (on.domain.end - on.domain.start);
      });
    this.parametersU = options.parametersU ?? located(curvesU[0]!, curvesV);
    this.parametersV = options.parametersV ?? located(curvesV[0]!, curvesU);
    if (this.parametersU.length !== curvesV.length || this.parametersV.length !== curvesU.length) {
      throw new RangeError('There must be one parameter per crossing curve.');
    }
    for (const params of [this.parametersU, this.parametersV]) {
      if (params[0] !== 0 || params[params.length - 1] !== 1) {
        throw new RangeError('Gordon parameters must run from 0 to 1.');
      }
    }
    this.blendU = new CardinalSplines(this.parametersU);
    this.blendV = new CardinalSplines(this.parametersV);
    this.families = { u, v };
    this.crossings = u.map((c, i) => this.parametersU.map((p, j) => {
      const point = c.at(p, 0);
      const other = v[j]!.at(this.parametersV[i]!, 0);
      const size = Math.max(curvesU[0]!.boundingBox().diagonal, curvesV[0]!.boundingBox().diagonal);
      if (point.distanceTo(other) > (options.tolerance ?? 1e-6 * size)) {
        throw new RangeError(`Curves ${i} and ${j} do not cross at their network parameters.`);
      }
      return point;
    }));
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.evaluateDerivative(u, v, 0, 0);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const weightsU = this.blendU.weights(u, orderU);
    const weightsV = this.blendV.weights(v, orderV);
    let sum = V3.ZERO;
    this.families.u.forEach((curve, i) => {
      sum = sum.addScaled(curve.at(u, orderU), weightsV[i]!);
    });
    this.families.v.forEach((curve, j) => {
      sum = sum.addScaled(curve.at(v, orderV), weightsU[j]!);
    });
    this.crossings.forEach((row, i) => {
      row.forEach((p, j) => {
        sum = sum.addScaled(p, -weightsV[i]! * weightsU[j]!);
      });
    });
    return sum;
  }

  override breakpointsU(): number[] {
    return mergedBreaks(this.families.u, this.parametersU);
  }

  override breakpointsV(): number[] {
    return mergedBreaks(this.families.v, this.parametersV);
  }

  override transformed(transform: Transform3): GordonSurface {
    return new GordonSurface(this.curvesU.map((c) => c.transformed(transform)), this.curvesV.map((c) => c.transformed(transform)), {
      parametersU: this.parametersU,
      parametersV: this.parametersV,
    });
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      curvesU: this.curvesU.map((c) => c.toJSON()),
      curvesV: this.curvesV.map((c) => c.toJSON()),
      parametersU: [...this.parametersU],
      parametersV: [...this.parametersV],
    };
  }
}
