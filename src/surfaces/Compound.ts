import type { Transform3, Vector3 } from 'mc-curves';
import { uniqueSorted } from '../math/numeric.js';
import type { ISurface, SurfaceDomain, SurfaceJSON } from './ISurface.js';
import { SurfaceBase } from './SurfaceBase.js';

const JOIN_TOLERANCE = 1e-7;

/**
 * Patches joined into a grid, indexed `[u][v]`. Each keeps its parameter speed, so the domain sizes
 * add up and derivatives equal those of the patches. The constructor throws on a gap.
 */
export class CompositeSurface extends SurfaceBase {
  readonly type = 'CompositeSurface';
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;
  private readonly startsU: number[];
  private readonly startsV: number[];

  constructor(readonly patches: readonly (readonly ISurface[])[]) {
    super();
    const countV = patches[0]?.length ?? 0;
    if (patches.length === 0 || countV === 0 || patches.some((column) => column.length !== countV)) {
      throw new RangeError('Composite patches must form a full grid.');
    }
    const widthU = (p: ISurface) => p.domainU.end - p.domainU.start;
    const widthV = (p: ISurface) => p.domainV.end - p.domainV.start;
    const relative = 1e-9;
    patches.forEach((column, i) => {
      column.forEach((patch, j) => {
        if (Math.abs(widthU(patch) - widthU(column[0]!)) > relative * widthU(column[0]!)) {
          throw new RangeError('Patches in one column must have the same u domain size.');
        }
        if (Math.abs(widthV(patch) - widthV(patches[0]![j]!)) > relative * widthV(patches[0]![j]!)) {
          throw new RangeError('Patches in one row must have the same v domain size.');
        }
        if (i > 0) {
          this.checkJoin(patches[i - 1]![j]!, patch, 'u');
        }
        if (j > 0) {
          this.checkJoin(column[j - 1]!, patch, 'v');
        }
      });
    });
    this.startsU = [0];
    for (const column of patches) {
      this.startsU.push(this.startsU[this.startsU.length - 1]! + widthU(column[0]!));
    }
    this.startsV = [0];
    for (const patch of patches[0]!) {
      this.startsV.push(this.startsV[this.startsV.length - 1]! + widthV(patch));
    }
    this.domainU = { start: 0, end: this.startsU[this.startsU.length - 1]! };
    this.domainV = { start: 0, end: this.startsV[this.startsV.length - 1]! };
  }

  private checkJoin(a: ISurface, b: ISurface, direction: 'u' | 'v'): void {
    const size = Math.max(a.boundingBox().diagonal, b.boundingBox().diagonal, 1);
    for (let k = 0; k <= 8; k++) {
      const s = k / 8;
      const [p, q] =
        direction === 'u'
          ? [a.pointAt(a.domainU.end, a.domainV.start + (a.domainV.end - a.domainV.start) * s), b.pointAt(b.domainU.start, b.domainV.start + (b.domainV.end - b.domainV.start) * s)]
          : [a.pointAt(a.domainU.start + (a.domainU.end - a.domainU.start) * s, a.domainV.end), b.pointAt(b.domainU.start + (b.domainU.end - b.domainU.start) * s, b.domainV.start)];
      if (p.distanceTo(q) > JOIN_TOLERANCE * size) {
        throw new RangeError('Neighbouring composite patches must share their edge.');
      }
    }
  }

  private locate(t: number, starts: readonly number[]): number {
    let lo = 0;
    let hi = starts.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (t < starts[mid]!) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return lo;
  }

  private local(u: number, v: number): { patch: ISurface; lu: number; lv: number } {
    const i = this.locate(u, this.startsU);
    const j = this.locate(v, this.startsV);
    const patch = this.patches[i]![j]!;
    return { patch, lu: patch.domainU.start + (u - this.startsU[i]!), lv: patch.domainV.start + (v - this.startsV[j]!) };
  }

  protected evaluate(u: number, v: number): Vector3 {
    const { patch, lu, lv } = this.local(u, v);
    return patch.pointAt(lu, lv);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    const { patch, lu, lv } = this.local(u, v);
    return patch.derivativeAt(lu, lv, orderU, orderV);
  }

  override breakpointsU(): number[] {
    const all = this.patches.flatMap((column, i) =>
      column.flatMap((patch) => patch.breakpointsU().map((b) => this.startsU[i]! + (b - patch.domainU.start))),
    );
    return uniqueSorted(all, 1e-10);
  }

  override breakpointsV(): number[] {
    const all = this.patches.flatMap((column) =>
      column.flatMap((patch, j) => patch.breakpointsV().map((b) => this.startsV[j]! + (b - patch.domainV.start))),
    );
    return uniqueSorted(all, 1e-10);
  }

  override transformed(transform: Transform3): CompositeSurface {
    return new CompositeSurface(this.patches.map((column) => column.map((patch) => patch.transformed(transform))));
  }

  toJSON(): SurfaceJSON {
    return { type: this.type, patches: this.patches.map((column) => column.map((patch) => patch.toJSON())) };
  }
}

export interface ParametricSurfaceOptions {
  /** Exact derivative function. Without it, derivatives are numeric. */
  readonly derivative?: (u: number, v: number, orderU: number, orderV: number) => Vector3;
  /** Interior parameter values where the surface may lose smoothness. */
  readonly breakpointsU?: readonly number[];
  readonly breakpointsV?: readonly number[];
  readonly periodicU?: boolean;
  readonly periodicV?: boolean;
}

/** Any function of two parameters. It holds a function, so it can't be saved to JSON. */
export class ParametricSurface extends SurfaceBase {
  readonly type = 'ParametricSurface';
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;

  constructor(
    readonly fn: (u: number, v: number) => Vector3,
    uStart: number,
    uEnd: number,
    vStart: number,
    vEnd: number,
    readonly options: ParametricSurfaceOptions = {},
  ) {
    super();
    if (!(uEnd > uStart) || !(vEnd > vStart)) {
      throw new RangeError('Surface domain must be of positive size.');
    }
    this.domainU = { start: uStart, end: uEnd };
    this.domainV = { start: vStart, end: vEnd };
  }

  override get isPeriodicU(): boolean {
    return this.options.periodicU ?? false;
  }

  override get isPeriodicV(): boolean {
    return this.options.periodicV ?? false;
  }

  protected evaluate(u: number, v: number): Vector3 {
    return this.fn(u, v);
  }

  protected override evaluateDerivative(u: number, v: number, orderU: number, orderV: number): Vector3 {
    return this.options.derivative ? this.options.derivative(u, v, orderU, orderV) : this.numericDerivative(u, v, orderU, orderV);
  }

  override breakpointsU(): number[] {
    return [this.domainU.start, ...uniqueSorted(this.options.breakpointsU ?? []).filter((b) => b > this.domainU.start && b < this.domainU.end), this.domainU.end];
  }

  override breakpointsV(): number[] {
    return [this.domainV.start, ...uniqueSorted(this.options.breakpointsV ?? []).filter((b) => b > this.domainV.start && b < this.domainV.end), this.domainV.end];
  }

  toJSON(): SurfaceJSON {
    throw new TypeError('A ParametricSurface holds a function, so it cannot be saved to JSON.');
  }
}
