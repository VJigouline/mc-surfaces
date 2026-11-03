import { Vector3 } from 'meercad-curves';

/** Distance below which two points count as coincident, at unit model size. */
export const GEOMETRIC_TOLERANCE = 1e-9;

export function clamp(v: number, lo: number, hi: number): number {
  return v < lo ? lo : v > hi ? hi : v;
}

/** Eight-point Gauss-Legendre rule on [-1, 1]: positive nodes and their weights. */
export const GL_NODES = [0.1834346424956498, 0.525532409916329, 0.7966664774136267, 0.9602898564975363];
export const GL_WEIGHTS = [0.362683783378362, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763];

/** Bisection root of `f` on [a, b]. The ends must have opposite signs, or one must be a root. */
export function bisect(f: (x: number) => number, a: number, b: number, iterations = 50): number {
  let lo = a;
  let hi = b;
  const fLo = f(lo);
  if (fLo === 0) {
    return lo;
  }
  for (let i = 0; i < iterations; i++) {
    const mid = (lo + hi) / 2;
    const fMid = f(mid);
    if (fMid === 0) {
      return mid;
    }
    if (fLo < 0 === fMid < 0) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return (lo + hi) / 2;
}

/** Golden-section maximum of a unimodal `f` on [a, b]. */
export function goldenMaximum(f: (x: number) => number, a: number, b: number, iterations = 30): number {
  const ratio = (Math.sqrt(5) - 1) / 2;
  let lo = a;
  let hi = b;
  let x1 = hi - ratio * (hi - lo);
  let x2 = lo + ratio * (hi - lo);
  let f1 = f(x1);
  let f2 = f(x2);
  for (let i = 0; i < iterations; i++) {
    if (f1 >= f2) {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = hi - ratio * (hi - lo);
      f1 = f(x1);
    } else {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + ratio * (hi - lo);
      f2 = f(x2);
    }
  }
  return (lo + hi) / 2;
}

export function binomial(n: number, k: number): number {
  let result = 1;
  for (let i = 1; i <= k; i++) {
    result = (result * (n - k + i)) / i;
  }
  return Math.round(result);
}

/**
 * Finite-difference weights for the `order`-th derivative at `x` over arbitrary `nodes`
 * (Fornberg's algorithm). Works for off-centre nodes, so it stays accurate at span ends.
 */
export function finiteDifferenceWeights(x: number, nodes: readonly number[], order: number): number[] {
  const n = nodes.length;
  const c: number[][] = Array.from({ length: n }, () => new Array<number>(order + 1).fill(0));
  let c1 = 1;
  let c4 = nodes[0]! - x;
  c[0]![0] = 1;
  for (let i = 1; i < n; i++) {
    const mn = Math.min(i, order);
    let c2 = 1;
    const c5 = c4;
    c4 = nodes[i]! - x;
    for (let j = 0; j < i; j++) {
      const c3 = nodes[i]! - nodes[j]!;
      c2 *= c3;
      if (j === i - 1) {
        for (let k = mn; k >= 1; k--) {
          c[i]![k] = (c1 * (k * c[i - 1]![k - 1]! - c5 * c[i - 1]![k]!)) / c2;
        }
        c[i]![0] = (-c1 * c5 * c[i - 1]![0]!) / c2;
      }
      for (let k = mn; k >= 1; k--) {
        c[j]![k] = (c4 * c[j]![k]! - k * c[j]![k - 1]!) / c3;
      }
      c[j]![0] = (c4 * c[j]![0]!) / c3;
    }
    c1 = c2;
  }
  return c.map((row) => row[order]!);
}

/** The `k`-th derivative of `cos t`. */
export function dCos(t: number, k: number): number {
  return k === 0 ? Math.cos(t) : Math.cos(t + (k * Math.PI) / 2);
}

/** The `k`-th derivative of `sin t`. */
export function dSin(t: number, k: number): number {
  return k === 0 ? Math.sin(t) : Math.sin(t + (k * Math.PI) / 2);
}

/** The `k`-th derivative of `cosh t`. */
export function dCosh(t: number, k: number): number {
  return k % 2 === 0 ? Math.cosh(t) : Math.sinh(t);
}

/** The `k`-th derivative of `sinh t`. */
export function dSinh(t: number, k: number): number {
  return k % 2 === 0 ? Math.sinh(t) : Math.cosh(t);
}

/** The `k`-th derivative of `a + b t`. */
export function dLinear(a: number, b: number, t: number, k: number): number {
  return k === 0 ? a + b * t : k === 1 ? b : 0;
}

/** Solves `A x = b` for number right-hand sides, by Gaussian elimination with partial pivoting. */
export function solveLinear(matrix: readonly (readonly number[])[], rhs: readonly number[]): number[] {
  const n = matrix.length;
  const a = matrix.map((row, i) => [...row, rhs[i]!]);
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let r = col + 1; r < n; r++) {
      if (Math.abs(a[r]![col]!) > Math.abs(a[pivot]![col]!)) {
        pivot = r;
      }
    }
    if (Math.abs(a[pivot]![col]!) < 1e-14) {
      throw new RangeError('Linear system is singular.');
    }
    [a[col], a[pivot]] = [a[pivot]!, a[col]!];
    for (let r = col + 1; r < n; r++) {
      const factor = a[r]![col]! / a[col]![col]!;
      if (factor !== 0) {
        for (let c = col; c <= n; c++) {
          a[r]![c]! -= factor * a[col]![c]!;
        }
      }
    }
  }
  const x = new Array<number>(n);
  for (let r = n - 1; r >= 0; r--) {
    let v = a[r]![n]!;
    for (let c = r + 1; c < n; c++) {
      v -= a[r]![c]! * x[c]!;
    }
    x[r] = v / a[r]![r]!;
  }
  return x;
}

/** Solves `A x = b` for vector right-hand sides. */
export function solveVectorSystem(matrix: readonly (readonly number[])[], rhs: readonly Vector3[]): Vector3[] {
  const xs = solveLinear(matrix, rhs.map((p) => p.x));
  const ys = solveLinear(matrix, rhs.map((p) => p.y));
  const zs = solveLinear(matrix, rhs.map((p) => p.z));
  return xs.map((x, i) => new Vector3(x, ys[i]!, zs[i]!));
}

/** Sorted values with near-duplicates removed. */
export function uniqueSorted(values: readonly number[], tolerance = 1e-12): number[] {
  const sorted = [...values].sort((a, b) => a - b);
  const result: number[] = [];
  for (const v of sorted) {
    if (result.length === 0 || v - result[result.length - 1]! > tolerance) {
      result.push(v);
    }
  }
  return result;
}

/** Points per smooth span, so that the whole domain gets at least `total` intervals. */
export function gridParameters(breaks: readonly number[], total: number): number[] {
  const spans = breaks.length - 1;
  const perSpan = spans > 100 ? 2 : Math.max(4, Math.ceil(total / spans));
  const result: number[] = [];
  for (let i = 0; i < spans; i++) {
    for (let k = 0; k < perSpan; k++) {
      result.push(breaks[i]! + ((breaks[i + 1]! - breaks[i]!) * k) / perSpan);
    }
  }
  result.push(breaks[spans]!);
  return result;
}
