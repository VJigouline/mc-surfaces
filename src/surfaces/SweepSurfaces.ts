import { Vector3, type ICurve } from 'meercad-curves';
import type { SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { SurfaceBase } from './SurfaceBase.js';

interface MovingFrame {
  readonly tangent: Vector3;
  /** Perpendicular to the tangent. */
  readonly normal: Vector3;
  /** Tangent cross normal. */
  readonly binormal: Vector3;
}

const TABLE_SIZE = 256;

/**
 * Frames along a curve that turn as little as possible about its tangent, by the double-reflection
 * method of Wang et al. (2008). A table is built once, and each query steps from its nearest entry.
 */
export class RotationMinimizingFrames {
  private readonly parameters: number[] = [];
  private readonly frames: MovingFrame[] = [];

  constructor(readonly path: ICurve) {
    const breaks = path.breakpoints();
    const perSpan = Math.max(2, Math.ceil(TABLE_SIZE / (breaks.length - 1)));
    for (let i = 0; i + 1 < breaks.length; i++) {
      for (let k = 0; k < perSpan; k++) {
        this.parameters.push(breaks[i]! + ((breaks[i + 1]! - breaks[i]!) * k) / perSpan);
      }
    }
    this.parameters.push(breaks[breaks.length - 1]!);
    const first = this.parameters[0]!;
    const tangent = path.tangentAt(first);
    const hint = path.normalAt(first);
    const normal = hint.sub(tangent.scale(hint.dot(tangent)));
    const unit = normal.length > 1e-9 ? normal.normalized() : tangent.anyPerpendicular();
    this.frames.push({ tangent, normal: unit, binormal: tangent.cross(unit).normalized() });
    for (let i = 1; i < this.parameters.length; i++) {
      this.frames.push(this.step(this.parameters[i - 1]!, this.frames[i - 1]!, this.parameters[i]!));
    }
  }

  private step(from: number, frame: MovingFrame, to: number): MovingFrame {
    const v1 = this.path.pointAt(to).sub(this.path.pointAt(from));
    const c1 = v1.lengthSquared;
    const tangent = this.path.tangentAt(to);
    let reflected = frame.normal;
    let reflectedTangent = frame.tangent;
    if (c1 > 0) {
      reflected = frame.normal.addScaled(v1, (-2 / c1) * v1.dot(frame.normal));
      reflectedTangent = frame.tangent.addScaled(v1, (-2 / c1) * v1.dot(frame.tangent));
    }
    const v2 = tangent.sub(reflectedTangent);
    const c2 = v2.lengthSquared;
    const twisted = c2 > 1e-30 ? reflected.addScaled(v2, (-2 / c2) * v2.dot(reflected)) : reflected;
    const normal = twisted.sub(tangent.scale(twisted.dot(tangent))).normalized();
    return { tangent, normal, binormal: tangent.cross(normal).normalized() };
  }

  frameAt(t: number): MovingFrame {
    const table = this.parameters;
    let lo = 0;
    let hi = table.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (t < table[mid]!) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    if (t >= table[hi]!) {
      return this.frames[hi]!;
    }
    return t === table[lo]! ? this.frames[lo]! : this.step(table[lo]!, this.frames[lo]!, t);
  }
}

export interface SweepOptions {
  /**
   * `rotation-minimizing` (default) keeps the profile at a fixed attitude to the path, with as little
   * twist as possible. `translation` only moves the profile, so it never turns.
   */
  readonly mode?: 'rotation-minimizing' | 'translation';
}

/**
 * A profile carried along a path. `u` is the profile's parameter and `v` the path's. The profile is
 * taken as positioned at the start of the path, so `pointAt(u, pathStart)` is the profile itself.
 */
export class SweepSurface extends SurfaceBase {
  readonly type = 'SweepSurface';
  readonly mode: 'rotation-minimizing' | 'translation';
  private readonly frames: RotationMinimizingFrames | undefined;
  private readonly startPoint: Vector3;

  constructor(
    readonly profile: ICurve,
    readonly path: ICurve,
    options: SweepOptions = {},
  ) {
    super();
    this.mode = options.mode ?? 'rotation-minimizing';
    this.frames = this.mode === 'translation' ? undefined : new RotationMinimizingFrames(path);
    this.startPoint = path.startPoint;
  }

  get domainU(): SurfaceDomain {
    return this.profile.domain;
  }

  get domainV(): SurfaceDomain {
    return this.path.domain;
  }

  override get isPeriodicU(): boolean {
    return this.profile.isPeriodic;
  }

  protected evaluate(u: number, v: number): Vector3 {
    const offset = this.profile.pointAt(u).sub(this.startPoint);
    const origin = this.path.pointAt(v);
    if (!this.frames) {
      return origin.add(offset);
    }
    const start = this.frames.frameAt(this.path.domain.start);
    const now = this.frames.frameAt(v);
    return origin
      .addScaled(now.tangent, offset.dot(start.tangent))
      .addScaled(now.normal, offset.dot(start.normal))
      .addScaled(now.binormal, offset.dot(start.binormal));
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    if (this.frames) {
      return this.numericDerivative(u, v, orderU, orderV);
    }
    if (orderU > 0 && orderV > 0) {
      return Vector3.ZERO;
    }
    return orderU > 0 ? this.profile.derivativeAt(u, orderU) : this.path.derivativeAt(v, orderV);
  }

  override breakpointsU(): number[] {
    return this.profile.breakpoints();
  }

  override breakpointsV(): number[] {
    return this.path.breakpoints();
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, profile: this.profile.toJSON(), path: this.path.toJSON(), mode: this.mode };
  }
}

/**
 * A tube of constant radius around a path: `P(v) + r (cos u N(v) + sin u B(v))`, with `u` the angle
 * in radians and `(N, B)` rotation-minimizing. It folds over itself where the radius exceeds the
 * path's radius of curvature.
 */
export class PipeSurface extends SurfaceBase {
  readonly type = 'PipeSurface';
  readonly domainU: SurfaceDomain = { start: 0, end: 2 * Math.PI };
  private readonly frames: RotationMinimizingFrames;

  constructor(
    readonly path: ICurve,
    readonly radius: number,
  ) {
    super();
    if (!(radius > 0) || !Number.isFinite(radius)) {
      throw new RangeError('Radius must be positive.');
    }
    this.frames = new RotationMinimizingFrames(path);
  }

  get domainV(): SurfaceDomain {
    return this.path.domain;
  }

  override get isPeriodicU(): boolean {
    return true;
  }

  protected evaluate(u: number, v: number): Vector3 {
    const f = this.frames.frameAt(v);
    return this.path.pointAt(v).addScaled(f.normal, this.radius * Math.cos(u)).addScaled(f.binormal, this.radius * Math.sin(u));
  }

  override breakpointsV(): number[] {
    return this.path.breakpoints();
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, path: this.path.toJSON(), radius: this.radius };
  }
}
