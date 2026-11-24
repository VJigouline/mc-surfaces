import { Vector3, type Transform3 } from 'mc-curves';
import type { SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { NurbsSurface } from './NurbsSurface.js';
import { SurfaceBase } from './SurfaceBase.js';

/**
 * Finds the cell holding `t`, right-handed at interior knots.
 * The last cell also holds the domain end.
 */
function cellOf(t: number, knots: readonly number[]): number {
  let lo = 0;
  let hi = knots.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (t < knots[mid]!) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return lo;
}

function checkGrid(points: readonly (readonly Vector3[])[], minimum: number): void {
  const countV = points[0]?.length ?? 0;
  if (points.length < minimum || countV < minimum) {
    throw new RangeError(`The grid needs at least ${minimum} by ${minimum} points.`);
  }
  if (points.some((row) => row.length !== countV)) {
    throw new RangeError('Every row of the grid must have the same length.');
  }
}

/**
 * Piecewise bilinear surface through a grid of points, indexed `[u][v]`: the surface counterpart of
 * a polyline. The parameters are the grid indices, so `pointAt(2.5, 0)` is halfway along an edge.
 */
export class GridSurface extends SurfaceBase {
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;
  readonly points: readonly (readonly Vector3[])[];

  constructor(points: readonly (readonly Vector3[])[]) {
    super();
    checkGrid(points, 2);
    this.points = points.map((row) => [...row]);
    this.domainU = { start: 0, end: points.length - 1 };
    this.domainV = { start: 0, end: points[0]!.length - 1 };
  }

  get type(): string {
    return 'GridSurface';
  }

  override get degreeU(): number {
    return 1;
  }

  override get degreeV(): number {
    return 1;
  }

  private cell(t: number, end: number): [number, number] {
    const i = Math.min(Math.floor(t), end - 1);
    return [i, t - i];
  }

  protected evaluate(u: number, v: number): Vector3 {
    const [i, a] = this.cell(u, this.domainU.end);
    const [j, b] = this.cell(v, this.domainV.end);
    const near = this.points[i]![j]!.lerp(this.points[i]![j + 1]!, b);
    const far = this.points[i + 1]![j]!.lerp(this.points[i + 1]![j + 1]!, b);
    return near.lerp(far, a);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    if (orderU > 1 || orderV > 1) {
      return Vector3.ZERO;
    }
    const [i, a] = this.cell(u, this.domainU.end);
    const [j, b] = this.cell(v, this.domainV.end);
    const [p00, p01, p10, p11] = [this.points[i]![j]!, this.points[i]![j + 1]!, this.points[i + 1]![j]!, this.points[i + 1]![j + 1]!];
    if (orderU === 1 && orderV === 1) {
      return p11.sub(p10).sub(p01).add(p00);
    }
    if (orderU === 1) {
      return p10.lerp(p11, b).sub(p00.lerp(p01, b));
    }
    return p01.lerp(p11, a).sub(p00.lerp(p10, a));
  }

  override breakpointsU(): number[] {
    return Array.from({ length: this.points.length }, (_, i) => i);
  }

  override breakpointsV(): number[] {
    return Array.from({ length: this.points[0]!.length }, (_, j) => j);
  }

  override transformed(transform: Transform3): GridSurface {
    return new GridSurface(this.points.map((row) => row.map((p) => transform.applyToPoint(p))));
  }

  override toNurbs(_tolerance?: number): NurbsSurface {
    const knots = (count: number) => [0, ...Array.from({ length: count }, (_, i) => i), count - 1];
    return new NurbsSurface(1, 1, this.points, knots(this.points.length), knots(this.points[0]!.length));
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, points: this.points.map((row) => row.map((p) => p.toArray())) };
  }
}

/** The patch through four corners, over [0, 1] in both directions. Twisted corners make it a saddle. */
export class BilinearSurface extends GridSurface {
  constructor(p00: Vector3, p10: Vector3, p01: Vector3, p11: Vector3) {
    super([
      [p00, p01],
      [p10, p11],
    ]);
  }

  override get type(): string {
    return 'BilinearSurface';
  }

  override transformed(transform: Transform3): BilinearSurface {
    const [[a, b], [c, d]] = this.points.map((row) => row.map((p) => transform.applyToPoint(p))) as [[Vector3, Vector3], [Vector3, Vector3]];
    return new BilinearSurface(a, c, b, d);
  }

  override toJSON(): SurfaceJSON {
    return {
      type: this.type,
      corners: [this.points[0]![0]!, this.points[1]![0]!, this.points[0]![1]!, this.points[1]![1]!].map((p) => p.toArray()),
    };
  }
}

/** The four cubic Hermite blending functions, or their `order`-th derivative with respect to `t`. */
function hermiteBasis(x: number, h: number, order: number): [number, number, number, number] {
  const scale = h ** -order;
  switch (order) {
    case 0:
      return [(2 * x - 3) * x * x + 1, h * ((x - 2) * x + 1) * x, (3 - 2 * x) * x * x, h * (x - 1) * x * x];
    case 1:
      return [scale * 6 * x * (x - 1), scale * h * (3 * x * x - 4 * x + 1), scale * 6 * x * (1 - x), scale * h * x * (3 * x - 2)];
    case 2:
      return [scale * (12 * x - 6), scale * h * (6 * x - 4), scale * (6 - 12 * x), scale * h * (6 * x - 2)];
    case 3:
      return [scale * 12, scale * h * 6, scale * -12, scale * h * 6];
    default:
      return [0, 0, 0, 0];
  }
}

export interface HermiteSurfaceData {
  /** Points, indexed `[u][v]`. */
  readonly points: readonly (readonly Vector3[])[];
  /** `Su` at each point. */
  readonly tangentsU: readonly (readonly Vector3[])[];
  /** `Sv` at each point. */
  readonly tangentsV: readonly (readonly Vector3[])[];
  /** `Suv` at each point. Default zero, which gives a Ferguson patch. */
  readonly twists?: readonly (readonly Vector3[])[];
  /** Default: 0, 1, 2, and so on. */
  readonly knotsU?: readonly number[];
  readonly knotsV?: readonly number[];
}

/** Bicubic patches from a point, two tangents and a twist at each grid point. C1 across the patches. */
export class HermiteSurface extends SurfaceBase {
  readonly type = 'HermiteSurface';
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;
  readonly points: readonly (readonly Vector3[])[];
  readonly tangentsU: readonly (readonly Vector3[])[];
  readonly tangentsV: readonly (readonly Vector3[])[];
  readonly twists: readonly (readonly Vector3[])[];
  readonly knotsU: readonly number[];
  readonly knotsV: readonly number[];

  constructor(data: HermiteSurfaceData) {
    super();
    const { points, tangentsU, tangentsV } = data;
    checkGrid(points, 2);
    const countU = points.length;
    const countV = points[0]!.length;
    const twists = data.twists ?? points.map((row) => row.map(() => Vector3.ZERO));
    for (const grid of [tangentsU, tangentsV, twists]) {
      if (grid.length !== countU || grid.some((row) => row.length !== countV)) {
        throw new RangeError('Tangents and twists must match the grid of points.');
      }
    }
    const knotsU = data.knotsU ?? Array.from({ length: countU }, (_, i) => i);
    const knotsV = data.knotsV ?? Array.from({ length: countV }, (_, j) => j);
    for (const [knots, count] of [[knotsU, countU], [knotsV, countV]] as const) {
      if (knots.length !== count) {
        throw new RangeError('There must be one knot per grid line.');
      }
      for (let i = 0; i + 1 < count; i++) {
        if (!(knots[i + 1]! > knots[i]!)) {
          throw new RangeError('Knots must be strictly increasing.');
        }
      }
    }
    this.points = points;
    this.tangentsU = tangentsU;
    this.tangentsV = tangentsV;
    this.twists = twists;
    this.knotsU = knotsU;
    this.knotsV = knotsV;
    this.domainU = { start: knotsU[0]!, end: knotsU[countU - 1]! };
    this.domainV = { start: knotsV[0]!, end: knotsV[countV - 1]! };
  }

  override get degreeU(): number {
    return 3;
  }

  override get degreeV(): number {
    return 3;
  }

  private blend(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const i = cellOf(u, this.knotsU);
    const j = cellOf(v, this.knotsV);
    const hu = this.knotsU[i + 1]! - this.knotsU[i]!;
    const hv = this.knotsV[j + 1]! - this.knotsV[j]!;
    const [a0, a1, a2, a3] = hermiteBasis((u - this.knotsU[i]!) / hu, hu, orderU);
    const [b0, b1, b2, b3] = hermiteBasis((v - this.knotsV[j]!) / hv, hv, orderV);
    let sum = Vector3.ZERO;
    for (const [ci, valueU, slopeU] of [[i, a0, a1], [i + 1, a2, a3]] as const) {
      for (const [cj, valueV, slopeV] of [[j, b0, b1], [j + 1, b2, b3]] as const) {
        sum = sum
          .addScaled(this.points[ci]![cj]!, valueU * valueV)
          .addScaled(this.tangentsU[ci]![cj]!, slopeU * valueV)
          .addScaled(this.tangentsV[ci]![cj]!, valueU * slopeV)
          .addScaled(this.twists[ci]![cj]!, slopeU * slopeV);
      }
    }
    return sum;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.blend(u, v, 0, 0);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.blend(u, v, orderU, orderV);
  }

  override breakpointsU(): number[] {
    return [...this.knotsU];
  }

  override breakpointsV(): number[] {
    return [...this.knotsV];
  }

  override transformed(transform: Transform3): HermiteSurface {
    const points = (grid: readonly (readonly Vector3[])[]) => grid.map((row) => row.map((p) => transform.applyToPoint(p)));
    const vectors = (grid: readonly (readonly Vector3[])[]) => grid.map((row) => row.map((d) => transform.applyToVector(d)));
    return new HermiteSurface({
      points: points(this.points),
      tangentsU: vectors(this.tangentsU),
      tangentsV: vectors(this.tangentsV),
      twists: vectors(this.twists),
      knotsU: this.knotsU,
      knotsV: this.knotsV,
    });
  }

  toJSON(): SurfaceJSON {
    const grid = (g: readonly (readonly Vector3[])[]) => g.map((row) => row.map((p) => p.toArray()));
    return {
      type: this.type,
      points: grid(this.points),
      tangentsU: grid(this.tangentsU),
      tangentsV: grid(this.tangentsV),
      twists: grid(this.twists),
      knotsU: [...this.knotsU],
      knotsV: [...this.knotsV],
    };
  }
}
