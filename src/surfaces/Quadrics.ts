import { BezierCurve, CircularArc, Frame3, HyperbolicArc, Vector3 } from 'mc-curves';
import { dCos, dCosh, dLinear, dSin, dSinh } from '../math/numeric.js';
import { FrameSurface, isFullTurn, type SurfaceRange } from './FrameSurface.js';
import { BezierSurface, revolveNurbs, type NurbsSurface } from './NurbsSurface.js';

const TURN: readonly [number, number] = [0, 2 * Math.PI];

function positive(value: number, name: string): void {
  if (!(value > 0) || !Number.isFinite(value)) {
    throw new RangeError(`${name} must be positive.`);
  }
}

/** A surface of revolution stretched afterwards: how the elliptic quadrics reach exact NURBS. */
function stretched(frame: Frame3, unit: NurbsSurface, a: number, b: number): NurbsSurface {
  return unit.mapControlPoints((p) => frame.toWorld(new Vector3(a * p.x, b * p.y, p.z)));
}

/** The half-plane through local X and Z, with X and Y of the frame along them. */
const XZ = Frame3.fromAxes(Vector3.ZERO, Vector3.X, Vector3.Z, Vector3.Y.negate());
/** The same half-plane, with the frame's X along local Z. */
const ZX = Frame3.fromAxes(Vector3.ZERO, Vector3.Z, Vector3.X, Vector3.Y);

/** `(a cos v cos u, b cos v sin u, c sin v)`. By default the whole ellipsoid. */
export class EllipsoidSurface extends FrameSurface {
  readonly type = 'EllipsoidSurface';

  constructor(
    frame: Frame3,
    readonly radiusX: number,
    readonly radiusY: number,
    readonly radiusZ: number,
    range: SurfaceRange = {},
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? [-Math.PI / 2, Math.PI / 2] });
    positive(radiusX, 'Radius X');
    positive(radiusY, 'Radius Y');
    positive(radiusZ, 'Radius Z');
    if (this.domainV.start < -Math.PI / 2 - 1e-12 || this.domainV.end > Math.PI / 2 + 1e-12) {
      throw new RangeError('Ellipsoid latitude must stay within a quarter turn of the equator.');
    }
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const ring = dCos(v, orderV);
    return new Vector3(
      this.radiusX * ring * dCos(u, orderU),
      this.radiusY * ring * dSin(u, orderU),
      orderU === 0 ? this.radiusZ * dSin(v, orderV) : 0,
    );
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): EllipsoidSurface {
    return new EllipsoidSurface(frame, this.radiusX * scale, this.radiusY * scale, this.radiusZ * scale, { u, v });
  }

  protected jsonFields(): Record<string, unknown> {
    return { radiusX: this.radiusX, radiusY: this.radiusY, radiusZ: this.radiusZ };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const profile = new CircularArc(XZ, 1, v0, v1).toNurbs();
    const unit = revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1);
    return unit.mapControlPoints((p) => this.frame.toWorld(new Vector3(this.radiusX * p.x, this.radiusY * p.y, this.radiusZ * p.z)));
  }
}

/** `(a v cos u, b v sin u, c v^2)`, with `v >= 0`. The vertex is at `v = 0`. */
export class EllipticParaboloidSurface extends FrameSurface {
  readonly type = 'EllipticParaboloidSurface';

  constructor(
    frame: Frame3,
    readonly a: number,
    readonly b: number,
    readonly c: number,
    range: SurfaceRange = {},
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? [0, 1] });
    positive(a, 'a');
    positive(b, 'b');
    positive(c, 'c');
    if (this.domainV.start < 0) {
      throw new RangeError('Paraboloid v must not be negative.');
    }
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const r = dLinear(0, 1, v, orderV);
    const height = orderU > 0 ? 0 : this.c * (orderV === 0 ? v * v : orderV === 1 ? 2 * v : orderV === 2 ? 2 : 0);
    return new Vector3(this.a * r * dCos(u, orderU), this.b * r * dSin(u, orderU), height);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): EllipticParaboloidSurface {
    return new EllipticParaboloidSurface(frame, this.a * scale, this.b * scale, this.c * scale, { u, v });
  }

  protected jsonFields(): Record<string, unknown> {
    return { a: this.a, b: this.b, c: this.c };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const profile = new BezierCurve([
      new Vector3(v0, 0, this.c * v0 * v0),
      new Vector3((v0 + v1) / 2, 0, this.c * v0 * v1),
      new Vector3(v1, 0, this.c * v1 * v1),
    ]).toNurbs();
    return stretched(this.frame, revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1), this.a, this.b);
  }
}

/** `scale (u, v, u^2 / a^2 - v^2 / b^2)`: the saddle. Its domain must be given. */
export class HyperbolicParaboloidSurface extends FrameSurface {
  readonly type = 'HyperbolicParaboloidSurface';

  constructor(
    frame: Frame3,
    readonly a: number,
    readonly b: number,
    range: Required<SurfaceRange>,
    readonly scale = 1,
  ) {
    super(frame, range);
    positive(a, 'a');
    positive(b, 'b');
    positive(scale, 'Scale');
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const s = this.scale;
    const square = (t: number, k: number) => (k === 0 ? t * t : k === 1 ? 2 * t : k === 2 ? 2 : 0);
    const x = orderV === 0 ? dLinear(0, s, u, orderU) : 0;
    const y = orderU === 0 ? dLinear(0, s, v, orderV) : 0;
    let z = 0;
    if (orderV === 0) {
      z += (s * square(u, orderU)) / this.a ** 2;
    }
    if (orderU === 0) {
      z -= (s * square(v, orderV)) / this.b ** 2;
    }
    return new Vector3(x, y, z);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): HyperbolicParaboloidSurface {
    return new HyperbolicParaboloidSurface(frame, this.a, this.b, { u, v }, this.scale * scale);
  }

  protected jsonFields(): Record<string, unknown> {
    return { a: this.a, b: this.b, scale: this.scale };
  }

  override get degreeU(): number {
    return 2;
  }

  override get degreeV(): number {
    return 2;
  }

  /** A biquadratic Bezier patch, with the same parameters once rescaled to [0, 1]. */
  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const s = this.scale;
    const xs = [u0, (u0 + u1) / 2, u1];
    const ys = [v0, (v0 + v1) / 2, v1];
    const f = [u0 * u0, u0 * u1, u1 * u1].map((w) => w / this.a ** 2);
    const g = [v0 * v0, v0 * v1, v1 * v1].map((w) => w / this.b ** 2);
    return new BezierSurface(
      xs.map((x, i) => ys.map((y, j) => this.frame.toWorld(new Vector3(s * x, s * y, s * (f[i]! - g[j]!))))),
    );
  }
}

/**
 * One sheet: `(a cosh v cos u, b cosh v sin u, c sinh v)`.
 * Two sheets, upper one: `(a sinh v cos u, b sinh v sin u, c cosh v)`, with `v >= 0` and its vertex at `v = 0`.
 */
export class HyperboloidSurface extends FrameSurface {
  readonly type = 'HyperboloidSurface';

  constructor(
    frame: Frame3,
    readonly a: number,
    readonly b: number,
    readonly c: number,
    range: SurfaceRange = {},
    readonly sheets: 1 | 2 = 1,
  ) {
    super(frame, { u: range.u ?? TURN, v: range.v ?? (sheets === 1 ? [-1, 1] : [0, 1]) });
    positive(a, 'a');
    positive(b, 'b');
    positive(c, 'c');
    if (sheets === 2 && this.domainV.start < 0) {
      throw new RangeError('On a two-sheet hyperboloid, v must not be negative.');
    }
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const one = this.sheets === 1;
    const ring = one ? dCosh(v, orderV) : dSinh(v, orderV);
    const height = orderU > 0 ? 0 : this.c * (one ? dSinh(v, orderV) : dCosh(v, orderV));
    return new Vector3(this.a * ring * dCos(u, orderU), this.b * ring * dSin(u, orderU), height);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): HyperboloidSurface {
    return new HyperboloidSurface(frame, this.a * scale, this.b * scale, this.c * scale, { u, v }, this.sheets);
  }

  protected jsonFields(): Record<string, unknown> {
    return { a: this.a, b: this.b, c: this.c, sheets: this.sheets };
  }

  override get isPeriodicU(): boolean {
    return isFullTurn(this.domainU);
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const [u0, u1] = this.rangeU;
    const [v0, v1] = this.rangeV;
    const profile =
      this.sheets === 1
        ? new HyperbolicArc(XZ, 1, this.c, v0, v1).toNurbs()
        : new HyperbolicArc(ZX, this.c, 1, v0, v1).toNurbs();
    return stretched(this.frame, revolveNurbs(profile, Vector3.ZERO, Vector3.Z, u0, u1), this.a, this.b);
  }
}

/**
 * `(scale v cos u, scale v sin u, pitch u / 2pi)`: a straight line turning about and climbing the axis,
 * as on a screw thread. `u` is the angle in radians and may span several turns. Its domain must be given.
 */
export class HelicoidSurface extends FrameSurface {
  readonly type = 'HelicoidSurface';

  constructor(
    frame: Frame3,
    readonly pitch: number,
    range: Required<SurfaceRange>,
    readonly scale = 1,
  ) {
    super(frame, range);
    if (!Number.isFinite(pitch) || pitch === 0) {
      throw new RangeError('Helicoid pitch must not be zero.');
    }
    positive(scale, 'Scale');
  }

  protected local(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const r = this.scale * dLinear(0, 1, v, orderV);
    const climb = orderV > 0 ? 0 : (this.pitch / (2 * Math.PI)) * dLinear(0, 1, u, orderU);
    return new Vector3(r * dCos(u, orderU), r * dSin(u, orderU), climb);
  }

  protected rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): HelicoidSurface {
    return new HelicoidSurface(frame, this.pitch * scale, { u, v }, this.scale * scale);
  }

  protected jsonFields(): Record<string, unknown> {
    return { pitch: this.pitch, scale: this.scale };
  }

  override get degreeV(): number {
    return 1;
  }
}
