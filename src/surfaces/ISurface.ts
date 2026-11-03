import type { BoundingBox3, ICurve, Transform3, Vector3 } from 'meercad-curves';
import type { NurbsSurface } from './NurbsSurface.js';

/** The closed interval of one parameter. */
export interface SurfaceDomain {
  readonly start: number;
  readonly end: number;
}

/** Parameters of a point on a surface. */
export interface SurfaceParameter {
  readonly u: number;
  readonly v: number;
}

/** Unit vectors at a surface point. `tangentV` is not generally perpendicular to `tangentU`. */
export interface SurfaceFrame {
  readonly origin: Vector3;
  readonly tangentU: Vector3;
  readonly tangentV: Vector3;
  readonly normal: Vector3;
}

/** First fundamental form: `E = Su.Su`, `F = Su.Sv`, `G = Sv.Sv`. */
export interface FirstFundamentalForm {
  readonly E: number;
  readonly F: number;
  readonly G: number;
}

/** Second fundamental form: `L = Suu.n`, `M = Suv.n`, `N = Svv.n`. */
export interface SecondFundamentalForm {
  readonly L: number;
  readonly M: number;
  readonly N: number;
}

/**
 * Principal curvatures, with `k1 >= k2`. A curvature is positive where the surface bends
 * towards its normal. The directions are unit tangents, and `direction2 = normal x direction1`.
 */
export interface PrincipalCurvatures {
  readonly k1: number;
  readonly k2: number;
  readonly direction1: Vector3;
  readonly direction2: Vector3;
}

export interface ClosestSurfacePoint extends SurfaceParameter {
  readonly point: Vector3;
  readonly distance: number;
}

export interface SurfaceTessellationOptions {
  /** Maximum distance between the surface and its triangles. Default 0.01. */
  readonly chordTolerance?: number;
  /** Maximum length of a grid edge. Default: unlimited. */
  readonly maxEdgeLength?: number;
  /** Refinement stops once the grid would exceed this many points. Default 1 000 000. */
  readonly maxPoints?: number;
}

/**
 * Triangle mesh. Points are not shared across a seam or at a pole, so each keeps its own
 * parameters, and triangle winding agrees with the surface normal.
 */
export interface SurfaceMesh {
  readonly points: readonly Vector3[];
  readonly normals: readonly Vector3[];
  readonly parameters: readonly SurfaceParameter[];
  /** Three point indices per triangle. */
  readonly triangles: readonly number[];
}

/** A crossing of a curve and a surface. */
export interface CurveSurfaceIntersection extends SurfaceParameter {
  readonly point: Vector3;
  /** Parameter on the curve; for `intersectLine`, the signed distance from the line origin. */
  readonly t: number;
  /** Distance between the curve point and the surface point. */
  readonly distance: number;
}

/** A polyline where a surface meets a plane. */
export interface SurfaceSection {
  readonly points: readonly Vector3[];
  readonly parameters: readonly SurfaceParameter[];
  /** True if the polyline forms a loop. The first point is not repeated at the end. */
  readonly closed: boolean;
}

/** A polyline where two surfaces meet. `parameters` is on this surface. */
export interface SurfaceIntersection extends SurfaceSection {
  readonly parametersOther: readonly SurfaceParameter[];
}

export interface SurfaceIntersectionOptions {
  /** Largest gap allowed at a result point. Default 1e-6. */
  readonly tolerance?: number;
  /** Chord tolerance of the mesh that seeds the search. Default: 1/300 of the bounding-box diagonal. */
  readonly meshTolerance?: number;
}

/** JSON form of a surface. Read it back with `surfaceFromJSON`. */
export interface SurfaceJSON {
  readonly type: string;
  readonly [key: string]: unknown;
}

/**
 * The interface every surface implements.
 *
 * Surfaces are immutable: operations such as `trim` and `transformed` return new surfaces. A surface
 * is a map from a finite rectangle of parameters to 3D space, so infinite planes are not represented.
 *
 * Parameters outside the domain are clamped, except in a periodic direction, where they wrap around.
 * The normal is always `Su x Sv`, so it depends on the parameterization. See `reversed`.
 */
export interface ISurface {
  /** Surface class name, such as `SphereSurface`. */
  readonly type: string;
  readonly domainU: SurfaceDomain;
  readonly domainV: SurfaceDomain;
  /** Polynomial degree in `u`, or `undefined` if the surface has none. */
  readonly degreeU: number | undefined;
  /** Polynomial degree in `v`, or `undefined` if the surface has none. */
  readonly degreeV: number | undefined;
  /** True if the edges at the two ends of the `u` domain coincide. */
  readonly isClosedU: boolean;
  /** True if the edges at the two ends of the `v` domain coincide. */
  readonly isClosedV: boolean;
  /** True if `u` wraps around its domain. */
  readonly isPeriodicU: boolean;
  /** True if `v` wraps around its domain. */
  readonly isPeriodicV: boolean;
  /** True if each direction is closed or has both end edges collapsed to points, as on a sphere. */
  readonly enclosesVolume: boolean;
  /** Total surface area. */
  readonly area: number;

  /** Point at the given parameters. */
  pointAt(u: number, v: number): Vector3;
  /**
   * Partial derivative with respect to `u` (`orderU` times) and `v` (`orderV` times).
   * Defaults to `Su`. At a joint it is the right-hand one in both parameters.
   */
  derivativeAt(u: number, v: number, orderU?: number, orderV?: number): Vector3;
  /** Unit `Su`, or a limit direction where `Su` vanishes. */
  tangentUAt(u: number, v: number): Vector3;
  /** Unit `Sv`, or a limit direction where `Sv` vanishes. */
  tangentVAt(u: number, v: number): Vector3;
  /** Unit normal `Su x Sv`. At a singular point such as a pole, its limit from inside the domain. */
  normalAt(u: number, v: number): Vector3;
  frameAt(u: number, v: number): SurfaceFrame;

  firstFundamentalFormAt(u: number, v: number): FirstFundamentalForm;
  secondFundamentalFormAt(u: number, v: number): SecondFundamentalForm;
  /** Principal curvatures. At a singular point, those at a nearby regular point. */
  principalCurvaturesAt(u: number, v: number): PrincipalCurvatures;
  /** Gaussian curvature: `k1 k2`. */
  gaussianCurvatureAt(u: number, v: number): number;
  /** Mean curvature: `(k1 + k2) / 2`. */
  meanCurvatureAt(u: number, v: number): number;

  /** `u` values where the surface may lose smoothness, including both domain ends. */
  breakpointsU(): number[];
  /** `v` values where the surface may lose smoothness, including both domain ends. */
  breakpointsV(): number[];
  /** The curve at fixed `u`, running in `v`. */
  isoCurveU(u: number): ICurve;
  /** The curve at fixed `v`, running in `u`. */
  isoCurveV(v: number): ICurve;
  /** The four edges: `v` start, `u` end, `v` end, `u` start. Each runs in its own increasing parameter. */
  boundaryCurves(): [ICurve, ICurve, ICurve, ICurve];

  /** Area of a parameter rectangle; the limits are clamped to the domain. */
  areaBetween(u0: number, u1: number, v0: number, v1: number): number;
  /** Volume enclosed if `enclosesVolume` is true, otherwise `undefined`. */
  enclosedVolume(): number | undefined;
  /** Centroid of the surface, treated as a thin sheet of uniform density. */
  centroid(): Vector3;

  /** A grid of `countU` by `countV` points, evenly spaced in the parameters, indexed `[u][v]`. */
  sample(countU: number, countV: number): Vector3[][];
  /** Triangle mesh within the chord tolerance. */
  tessellate(options?: SurfaceTessellationOptions): SurfaceMesh;
  boundingBox(): BoundingBox3;
  closestPoint(point: Vector3): ClosestSurfacePoint;
  distanceTo(point: Vector3): number;

  /** Where a curve crosses the surface. See the documentation for tangent contacts. */
  intersectCurve(curve: ICurve, tolerance?: number): CurveSurfaceIntersection[];
  /** Where an infinite line crosses the surface, sorted along `direction`. */
  intersectLine(origin: Vector3, direction: Vector3, tolerance?: number): CurveSurfaceIntersection[];
  /** Section by a plane, as polylines. A surface lying in the plane gives no section. */
  intersectPlane(origin: Vector3, normal: Vector3, options?: SurfaceIntersectionOptions): SurfaceSection[];
  /** Intersection polylines with another surface. See the documentation for overlaps and touches. */
  intersect(other: ISurface, options?: SurfaceIntersectionOptions): SurfaceIntersection[];

  /** True if all points lie in one plane within the tolerance. */
  isPlanar(tolerance?: number): boolean;
  /** Unit normal of the surface's plane, on the side of `normalAt`; `undefined` if not planar. */
  planeNormal(tolerance?: number): Vector3 | undefined;

  /** The piece over a parameter rectangle. It keeps this surface's parameters. */
  trim(u0: number, u1: number, v0: number, v1: number): ISurface;
  /** Splits at `u` into the pieces before and after. */
  splitU(u: number): [ISurface, ISurface];
  /** Splits at `v` into the pieces before and after. */
  splitV(v: number): [ISurface, ISurface];
  /** The same geometry with the normal flipped, by running `u` the other way over the same domain. */
  reversed(): ISurface;
  /** The surface under an affine transform. It keeps this surface's parameters. */
  transformed(transform: Transform3): ISurface;
  /** Parallel surface at signed `distance` along the normal. */
  offset(distance: number): ISurface;
  /**
   * NURBS form. Exact where the class supports it; otherwise a bicubic approximation
   * within `tolerance` (default 1e-4).
   */
  toNurbs(tolerance?: number): NurbsSurface;

  /** True if both surfaces agree at evenly spaced normalized parameters. */
  equals(other: ISurface, tolerance?: number): boolean;
  toJSON(): SurfaceJSON;
}
