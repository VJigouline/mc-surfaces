/**
 * B-spline basis functions. Algorithm numbers refer to Piegl and Tiller, The NURBS Book.
 */

/**
 * The knot span holding `u`, so `U[span] <= u < U[span + 1]`. At the domain end it is the last
 * non-empty span, which makes evaluation right-handed at interior knots.
 */
export function findSpan(u: number, degree: number, knots: readonly number[], count: number): number {
  if (u >= knots[count]!) {
    let span = count - 1;
    while (span > degree && knots[span]! >= knots[span + 1]!) {
      span--;
    }
    return span;
  }
  let lo = degree;
  let hi = count;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (u < knots[mid]!) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return lo;
}

/**
 * Non-zero basis functions at `u` and their derivatives (A2.3).
 * `result[k][j]` is the `k`-th derivative of the basis function `span - degree + j`.
 */
export function basisDerivatives(u: number, span: number, degree: number, knots: readonly number[], order: number): number[][] {
  const p = degree;
  const ndu: number[][] = Array.from({ length: p + 1 }, () => new Array<number>(p + 1).fill(0));
  const left = new Array<number>(p + 1).fill(0);
  const right = new Array<number>(p + 1).fill(0);
  ndu[0]![0] = 1;
  for (let j = 1; j <= p; j++) {
    left[j] = u - knots[span + 1 - j]!;
    right[j] = knots[span + j]! - u;
    let saved = 0;
    for (let r = 0; r < j; r++) {
      ndu[j]![r] = right[r + 1]! + left[j - r]!;
      const temp = ndu[r]![j - 1]! / ndu[j]![r]!;
      ndu[r]![j] = saved + right[r + 1]! * temp;
      saved = left[j - r]! * temp;
    }
    ndu[j]![j] = saved;
  }
  const result: number[][] = Array.from({ length: order + 1 }, () => new Array<number>(p + 1).fill(0));
  for (let j = 0; j <= p; j++) {
    result[0]![j] = ndu[j]![p]!;
  }
  const top = Math.min(order, p);
  for (let r = 0; r <= p; r++) {
    let s1 = 0;
    let s2 = 1;
    const a: number[][] = [new Array<number>(p + 1).fill(0), new Array<number>(p + 1).fill(0)];
    a[0]![0] = 1;
    for (let k = 1; k <= top; k++) {
      let d = 0;
      const rk = r - k;
      const pk = p - k;
      if (r >= k) {
        a[s2]![0] = a[s1]![0]! / ndu[pk + 1]![rk]!;
        d = a[s2]![0]! * ndu[rk]![pk]!;
      }
      const j1 = rk >= -1 ? 1 : -rk;
      const j2 = r - 1 <= pk ? k - 1 : p - r;
      for (let j = j1; j <= j2; j++) {
        a[s2]![j] = (a[s1]![j]! - a[s1]![j - 1]!) / ndu[pk + 1]![rk + j]!;
        d += a[s2]![j]! * ndu[rk + j]![pk]!;
      }
      if (r <= pk) {
        a[s2]![k] = -a[s1]![k - 1]! / ndu[pk + 1]![r]!;
        d += a[s2]![k]! * ndu[r]![pk]!;
      }
      result[k]![r] = d;
      [s1, s2] = [s2, s1];
    }
  }
  let factor = p;
  for (let k = 1; k <= top; k++) {
    for (let j = 0; j <= p; j++) {
      result[k]![j]! *= factor;
    }
    factor *= p - k;
  }
  return result;
}

/** Clamped knots with even interior spacing over [0, 1]. */
export function uniformClampedKnots(degree: number, count: number): number[] {
  const spans = count - degree;
  return [
    ...new Array<number>(degree).fill(0),
    ...Array.from({ length: spans + 1 }, (_, i) => i / spans),
    ...new Array<number>(degree).fill(1),
  ];
}

/** Knots by averaging the interpolation parameters (The NURBS Book, eq. 9.8). */
export function averagedKnots(params: readonly number[], degree: number): number[] {
  const n = params.length - 1;
  const knots = new Array<number>(n + degree + 2).fill(0);
  for (let j = 1; j <= n - degree; j++) {
    let sum = 0;
    for (let i = j; i < j + degree; i++) {
      sum += params[i]!;
    }
    knots[j + degree] = sum / degree;
  }
  for (let i = n + 1; i <= n + degree + 1; i++) {
    knots[i] = 1;
  }
  return knots;
}

/** One row of the interpolation matrix: all `count` basis functions at `u`. */
export function basisRow(u: number, degree: number, knots: readonly number[], count: number): number[] {
  const span = findSpan(u, degree, knots, count);
  const values = basisDerivatives(u, span, degree, knots, 0)[0]!;
  const row = new Array<number>(count).fill(0);
  values.forEach((value, j) => {
    row[span - degree + j] = value;
  });
  return row;
}
