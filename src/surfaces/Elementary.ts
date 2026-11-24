import { CircularArc, Frame3, Line, Vector3, type ICurve } from 'mc-curves';
import { clamp, dCos, dLinear, dSin } from '../math/numeric.js';
import { FrameSurface, isFullTurn, type SurfaceRange } from './FrameSurface.js';
import type { ClosestSurfacePoint, ISurface } from './ISurface.js';
import { BSplineSurface, revolveNurbs, type NurbsSurface } from './NurbsSurface.js';
import { OffsetSurface } from './SurfaceBase.js';

const TURN: readonly [number, number] = [0, 2 * Math.PI];

function positive(value: number, name: string): void {
  if (!(value > 0) || !Number.isFinite(value)) {
    throw new RangeError(`${name} must be positive.`);
  }
}

/** A frame whose X and Y axes are local X and Z: the half-plane a revolved profile lies in. */
const PROFILE_PLANE = Frame3.fromAxes(Vector3.ZERO, Vector3.X, Vector3.Z, Vector3.Y.negate());

/** A rectangular piece of a plane: `origin + scale (u X + v Y)`. */
export class PlaneSurface extends FrameSurface {
  readonly type = 'PlaneSurface';

  /** `scale` multiplies the parameters. It is how a scaled copy keeps its parameters. */
  constructor(
    frame: Frame3,
    range: Required<SurfaceRange>,
    readonly scale = 1,
  ) {
    super(frame, range);
    positive(scale, 'Scale');
  }

  /** A plane patch centred on `origin`, with the given sizes. */
  static centred(origin: Vector3, normal: Vector3, sizeU: number, sizeV: number, xHint?: Vector3): PlaneSurface {
    return new PlaneSurface(Frame3.fromNormal(origin, normal, xHint), { u: [-sizeU / 2, sizeU / 2], v: [-sizeV / 2, sizeV / 2] });
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    if (orderU === 0 && orderV === 0) {
      return new Vector3(this.scale * u, this.scale * v, 0);
    }
    if (orderU === 1 && orderV === 0) {
      return new Vector3(this.scale, 0, 0);
    }
    return orderU === 0 && orderV === 1 ? new Vector3(0, this.scale, 0) : Vector3.ZERO;
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): PlaneSurface {
    return new PlaneSurface(frame, { u, v }, this.scale * scale);
  }

  protected jsonFields(): Record<string, unknown> {
    return { scale: this.scale };
  }

  override get degreeU(): number {
    return 1;
  }

  override get degreeV(): number {
    return 1;
  }

  override normalAt(_u: number, _v: number): Vector3 {
    return this.frame.planeNormal;
  }

  override get area(): number {
    return (this.domainU.end - this.domainU.start) * (this.domainV.end - this.domainV.start) * this.scale ** 2;
  }

  override closestPoint(point: Vector3): ClosestSurfacePoint {
    const p = this.frame.toLocal(point);
    const u = clamp(p.x / this.scale, this.domainU.start, this.domainU.end);
    const v = clamp(p.y / this.scale, this.domainV.start, this.domainV.end);
    const q = this.pointAt(u, v);
    return { u, v, point: q, distance: q.distanceTo(point) };
  }

  override planeNormal(_tolerance?: number): Vector3 {
    return this.frame.planeNormal;
  }

  override isoCurveU(u: number): ICurve {
    return new Line(this.pointAt(u, this.domainV.start), this.pointAt(u, this.domainV.end));
  }

  override isoCurveV(v: number): ICurve {
    return new Line(this.pointAt(this.domainU.start, v), this.pointAt(this.domainU.end, v));
  }

  override offset(distance: number): ISurface {
    const f = this.frame;
    const moved = Frame3.fromAxes(f.origin.addScaled(f.planeNormal, distance), f.xAxis, f.yAxis, f.zAxis);
    return new PlaneSurface(moved, { u: this.rangeU, v: this.rangeV }, this.scale);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    return new BSplineSurface(
      1,
      1,
      [
        [this.pointAt(u0, v0), this.pointAt(u0, v1)],
        [this.pointAt(u1, v0), this.pointAt(u1, v1)],
      ],
      [u0, u0, u1, u1],
      [v0, v0, v1, v1],
    );
  }
}

/** `(r cos u, r sin u, scale v)`. By default a full turn, with `v` in [0, 1]. */
export class CylinderSurface extends FrameSurface {
  readonly type = 'CylinderSurface';

  constructor(
    frame: Frame3,
    readonly radius: number,
    range: SurfaceRange = {},
    readonly scale = 1,
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? [0, 1] });
    positive(radius, 'Radius');
    positive(scale, 'Scale');
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const r = orderV === 0 ? this.radius : 0;
    return new Vector3(r * dCos(u, orderU), r * dSin(u, orderU), orderU === 0 ? this.scale * dLinear(0, 1, v, orderV) : 0);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): CylinderSurface {
    return new CylinderSurface(frame, this.radius * scale, { u, v }, this.scale * scale);
  }

  protected jsonFields(): Record<string, unknown> {
    return { radius: this.radius, scale: this.scale };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override get area(): number {
    return (this.domainU.end - this.domainU.start) * (this.domainV.end - this.domainV.start) * this.radius * this.scale;
  }

  override closestPoint(point: Vector3): ClosestSurfacePoint {
    const p = this.frame.toLocal(point);
    const v = clamp(p.z / this.scale, this.domainV.start, this.domainV.end);
    let u = Math.atan2(p.y, p.x);
    const { start, end } = this.domainU;
    const span = 2 * Math.PI;
    u = start + ((((u - start) % span) + span) % span);
    if (u > end) {
      // Outside the arc, so the nearer end. On a point on the axis, every angle ties.
      const toEnd = u - end;
      const toStart = start + span - u;
      u = toEnd <= toStart ? end : start;
    }
    const q = this.pointAt(u, v);
    return { u, v, point: q, distance: q.distanceTo(point) };
  }

  override offset(distance: number): ISurface {
    const radius = this.radius + this.handedness * distance;
    if (!(radius > 0)) {
      return new OffsetSurface(this, distance);
    }
    return new CylinderSurface(this.frame, radius, { u: this.rangeU, v: this.rangeV }, this.scale);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const profile = new Line(new Vector3(this.radius, 0, this.scale * v0), new Vector3(this.radius, 0, this.scale * v1)).toNurbs();
    return revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1).mapControlPoints((p) => this.frame.toWorld(p));
  }
}

/**
 * `(r(v) cos u, r(v) sin u, scale v)`, with `r(v) = radius + scale v tan(semiAngle)`.
 * The radius must not go negative. It can reach zero, which is the apex.
 */
export class ConeSurface extends FrameSurface {
  readonly type = 'ConeSurface';

  constructor(
    frame: Frame3,
    readonly radius: number,
    readonly semiAngle: number,
    range: SurfaceRange = {},
    readonly scale = 1,
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? [0, 1] });
    positive(scale, 'Scale');
    if (!(Math.abs(semiAngle) < Math.PI / 2)) {
      throw new RangeError('Cone semi-angle must be below a right angle.');
    }
    const slope = this.slope;
    if (radius + slope * this.domainV.start < -1e-12 || radius + slope * this.domainV.end < -1e-12) {
      throw new RangeError('The cone radius must not go negative over the domain.');
    }
    if (!(Math.max(radius + slope * this.domainV.start, radius + slope * this.domainV.end) > 0)) {
      throw new RangeError('The cone radius must not be zero throughout.');
    }
  }

  private get slope(): number {
    return this.scale * Math.tan(this.semiAngle);
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const r = dLinear(this.radius, this.slope, v, orderV);
    return new Vector3(r * dCos(u, orderU), r * dSin(u, orderU), orderU === 0 ? this.scale * dLinear(0, 1, v, orderV) : 0);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): ConeSurface {
    return new ConeSurface(frame, this.radius * scale, this.semiAngle, { u, v }, this.scale * scale);
  }

  protected jsonFields(): Record<string, unknown> {
    return { radius: this.radius, semiAngle: this.semiAngle, scale: this.scale };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override get area(): number {
    const [v0, v1] = this.rangeV;
    const r0 = this.radius + this.slope * v0;
    const r1 = this.radius + this.slope * v1;
    return ((this.domainU.end - this.domainU.start) * (r0 + r1) * 0.5 * (this.scale * (v1 - v0))) / Math.cos(this.semiAngle);
  }

  /** Another cone with the same semi-angle. Its `v` shifts, because the normal leans along the axis. */
  override offset(distance: number): ISurface {
    const d = this.handedness * distance;
    const radius = this.radius + d / Math.cos(this.semiAngle);
    const shift = (d * Math.sin(this.semiAngle)) / this.scale;
    const v: [number, number] = [this.domainV.start - shift, this.domainV.end - shift];
    const slope = this.slope;
    if (radius + slope * v[0] < 0 || radius + slope * v[1] < 0) {
      return new OffsetSurface(this, distance);
    }
    return new ConeSurface(this.frame, radius, this.semiAngle, { u: this.rangeU, v }, this.scale);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const at = (v: number) => new Vector3(this.radius + this.slope * v, 0, this.scale * v);
    const profile = new Line(at(v0), at(v1)).toNurbs();
    return revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1).mapControlPoints((p) => this.frame.toWorld(p));
  }
}

/** `r (cos v cos u, cos v sin u, sin v)`: `u` is longitude and `v` latitude. By default the whole sphere. */
export class SphereSurface extends FrameSurface {
  readonly type = 'SphereSurface';

  constructor(
    frame: Frame3,
    readonly radius: number,
    range: SurfaceRange = {},
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? [-Math.PI / 2, Math.PI / 2] });
    positive(radius, 'Radius');
    if (this.domainV.start < -Math.PI / 2 - 1e-12 || this.domainV.end > Math.PI / 2 + 1e-12) {
      throw new RangeError('Sphere latitude must stay within a quarter turn of the equator.');
    }
  }

  /** A whole sphere with its frame aligned to the world axes. */
  static centredAt(centre: Vector3, radius: number): SphereSurface {
    return new SphereSurface(Frame3.fromNormal(centre), radius);
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const r = this.radius;
    const ring = r * dCos(v, orderV);
    return new Vector3(ring * dCos(u, orderU), ring * dSin(u, orderU), orderU === 0 ? r * dSin(v, orderV) : 0);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): SphereSurface {
    return new SphereSurface(frame, this.radius * scale, { u, v });
  }

  protected jsonFields(): Record<string, unknown> {
    return { radius: this.radius };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override get area(): number {
    const { domainU: du, domainV: dv } = this;
    return this.radius ** 2 * (du.end - du.start) * (Math.sin(dv.end) - Math.sin(dv.start));
  }

  override closestPoint(point: Vector3): ClosestSurfacePoint {
    const p = this.frame.toLocal(point);
    const longitude = Math.atan2(p.y, p.x);
    const latitude = Math.atan2(p.z, Math.hypot(p.x, p.y));
    const { start, end } = this.domainU;
    const u = start + ((((longitude - start) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI));
    if (u <= end && latitude >= this.domainV.start && latitude <= this.domainV.end) {
      const q = this.pointAt(u, latitude);
      return { u, v: latitude, point: q, distance: q.distanceTo(point) };
    }
    // Outside the patch, the answer lies on an edge, which needs the general search.
    return super.closestPoint(point);
  }

  override offset(distance: number): ISurface {
    const radius = this.radius + this.handedness * distance;
    if (!(radius > 0)) {
      return new OffsetSurface(this, distance);
    }
    return new SphereSurface(this.frame, radius, { u: this.rangeU, v: this.rangeV });
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const profile = new CircularArc(PROFILE_PLANE, this.radius, v0, v1).toNurbs();
    return revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1).mapControlPoints((p) => this.frame.toWorld(p));
  }
}

/** `((R + r cos v) cos u, (R + r cos v) sin u, r sin v)`. By default the whole torus. */
export class TorusSurface extends FrameSurface {
  readonly type = 'TorusSurface';

  constructor(
    frame: Frame3,
    readonly majorRadius: number,
    readonly minorRadius: number,
    range: SurfaceRange = {},
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? TURN });
    positive(majorRadius, 'Major radius');
    positive(minorRadius, 'Minor radius');
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const ring = orderV === 0 ? this.majorRadius + this.minorRadius * Math.cos(v) : this.minorRadius * dCos(v, orderV);
    return new Vector3(ring * dCos(u, orderU), ring * dSin(u, orderU), orderU === 0 ? this.minorRadius * dSin(v, orderV) : 0);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): TorusSurface {
    return new TorusSurface(frame, this.majorRadius * scale, this.minorRadius * scale, { u, v });
  }

  protected jsonFields(): Record<string, unknown> {
    return { majorRadius: this.majorRadius, minorRadius: this.minorRadius };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override get isPeriodicV(): boolean {
    return isFullTurn(this.domainV);
  }

  override offset(distance: number): ISurface {
    const minor = this.minorRadius + this.handedness * distance;
    if (!(minor > 0)) {
      return new OffsetSurface(this, distance);
    }
    return new TorusSurface(this.frame, this.majorRadius, minor, { u: this.rangeU, v: this.rangeV });
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const tube = Frame3.fromAxes(new Vector3(this.majorRadius, 0, 0), Vector3.X, Vector3.Z, Vector3.Y.negate());
    const profile = new CircularArc(tube, this.minorRadius, v0, v1).toNurbs();
    return revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1).mapControlPoints((p) => this.frame.toWorld(p));
  }
}
