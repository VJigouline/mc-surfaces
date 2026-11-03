import { Frame3, type Transform3, type Vector3 } from 'meercad-curves';
import type { ISurface, SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { SurfaceBase, TransformedSurface } from './SurfaceBase.js';

/** Parameter limits of a surface. A missing limit takes the class's default. */
export interface SurfaceRange {
  readonly u?: readonly [number, number];
  readonly v?: readonly [number, number];
}

const FULL_TURN = 2 * Math.PI;

/** True if an angle range covers exactly one turn. */
export function isFullTurn(domain: SurfaceDomain): boolean {
  return Math.abs(domain.end - domain.start - FULL_TURN) <= 1e-12 * FULL_TURN;
}

/**
 * Surface defined in closed form in a local coordinate system. Similarity transforms stay exact and
 * keep the parameters, by moving the frame and scaling the sizes; other affine transforms fall back
 * to `TransformedSurface`.
 */
export abstract class FrameSurface extends SurfaceBase {
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;

  constructor(
    readonly frame: Frame3,
    range: { u: readonly [number, number]; v: readonly [number, number] },
  ) {
    super();
    for (const [start, end] of [range.u, range.v]) {
      if (!Number.isFinite(start) || !Number.isFinite(end) || !(end! > start!)) {
        throw new RangeError('Surface domain must be finite and of positive size.');
      }
    }
    this.domainU = { start: range.u[0], end: range.u[1] };
    this.domainV = { start: range.v[0], end: range.v[1] };
  }

  /** Point (orders 0, 0) or derivative in local coordinates. */
  protected abstract local(u: number, v: number, orderU: number, orderV: number): Vector3;

  /** The same kind of surface in a new frame, with sizes multiplied by `scale`, over the given range. */
  protected abstract rebuild(frame: Frame3, scale: number, u: readonly [number, number], v: readonly [number, number]): FrameSurface;

  /** Surface-specific members of the JSON form. */
  protected abstract jsonFields(): Record<string, unknown>;

  protected get rangeU(): [number, number] {
    return [this.domainU.start, this.domainU.end];
  }

  protected get rangeV(): [number, number] {
    return [this.domainV.start, this.domainV.end];
  }

  /** +1 in a right-handed frame, -1 in a mirrored one. In a mirrored frame the normal turns inside out. */
  protected get handedness(): number {
    return this.frame.planeNormal.dot(this.frame.zAxis) < 0 ? -1 : 1;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.frame.toWorld(this.local(u, v, 0, 0));
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.frame.vectorToWorld(this.local(u, v, orderU, orderV));
  }

  override trim(u0: number, u1: number, v0: number, v1: number): ISurface {
    const [a, b, c, d] = this.checkedRange(u0, u1, v0, v1);
    return this.rebuild(this.frame, 1, [a, b], [c, d]);
  }

  override transformed(transform: Transform3): ISurface {
    const scale = transform.similarityScale();
    if (scale === undefined) {
      return new TransformedSurface(this, transform);
    }
    return this.rebuild(this.frame.transformed(transform), scale, this.rangeU, this.rangeV);
  }

  toJSON(): SurfaceJSON {
    return {
      type: this.type,
      frame: this.frame.toJSON(),
      u: this.rangeU,
      v: this.rangeV,
      ...this.jsonFields(),
    };
  }
}
