import { solveLinear } from './numeric.js';

/**
 * Cardinal cubic splines: function `j` is 1 at knot `j` and 0 at the others, and each is C2.
 * Any data through the knots is the sum of its values times these functions, which is how lofts
 * blend whole curves. Ends are natural, or periodic if requested. Two knots give straight lines.
 */
export class CardinalSplines {
  /** Number of independent nodes: the closing knot of a periodic set repeats the first. */
  readonly count: number;
  private readonly moments: number[][];

  constructor(
    readonly knots: readonly number[],
    readonly periodic = false,
  ) {
    const n = knots.length - 1;
    if (n < 1) {
      throw new RangeError('Cardinal splines need at least two knots.');
    }
    for (let i = 0; i < n; i++) {
      if (!(knots[i + 1]! > knots[i]!)) {
        throw new RangeError('Knots must be strictly increasing.');
      }
    }
    if (periodic && n < 3) {
      throw new RangeError('A closed set needs at least three different nodes.');
    }
    this.count = periodic ? n : n + 1;
    this.moments = Array.from({ length: this.count }, (_, j) => this.solveMoments(j));
  }

  /** Second derivatives at the knots of the function that is 1 at node `node`. */
  private solveMoments(node: number): number[] {
    const n = this.knots.length - 1;
    const h = (i: number) => this.knots[i + 1]! - this.knots[i]!;
    const y = (i: number) => (i % this.count === node ? 1 : 0);
    if (!this.periodic) {
      if (n === 1) {
        return [0, 0];
      }
      const size = n + 1;
      const a = Array.from({ length: size }, () => new Array<number>(size).fill(0));
      const b = new Array<number>(size).fill(0);
      a[0]![0] = 1;
      a[n]![n] = 1;
      for (let i = 1; i < n; i++) {
        a[i]![i - 1] = h(i - 1);
        a[i]![i] = 2 * (h(i - 1) + h(i));
        a[i]![i + 1] = h(i);
        b[i] = 6 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i - 1));
      }
      return solveLinear(a, b);
    }
    const a = Array.from({ length: n }, () => new Array<number>(n).fill(0));
    const b = new Array<number>(n).fill(0);
    for (let i = 0; i < n; i++) {
      const before = (i + n - 1) % n;
      const after = (i + 1) % n;
      a[i]![before]! += h(before);
      a[i]![i]! += 2 * (h(before) + h(i));
      a[i]![after]! += h(i);
      b[i] = 6 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(before)) / h(before));
    }
    const m = solveLinear(a, b);
    return [...m, m[0]!];
  }

  /** The `order`-th derivative of every cardinal function at `t`. Right-handed at knots. */
  weights(t: number, order = 0): number[] {
    const n = this.knots.length - 1;
    let i = 0;
    while (i < n - 1 && t >= this.knots[i + 1]!) {
      i++;
    }
    const h = this.knots[i + 1]! - this.knots[i]!;
    const x = t - this.knots[i]!;
    return this.moments.map((m, j) => {
      const y0 = i % this.count === j ? 1 : 0;
      const y1 = (i + 1) % this.count === j ? 1 : 0;
      const m0 = m[i]!;
      const m1 = m[i + 1]!;
      const b = (y1 - y0) / h - (h * (2 * m0 + m1)) / 6;
      const c = m0 / 2;
      const d = (m1 - m0) / (6 * h);
      switch (order) {
        case 0:
          return y0 + x * (b + x * (c + x * d));
        case 1:
          return b + x * (2 * c + 3 * x * d);
        case 2:
          return 2 * c + 6 * x * d;
        case 3:
          return 6 * d;
        default:
          return 0;
      }
    });
  }
}
