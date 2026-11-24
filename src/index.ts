// The shared math types come from the curve library, so one set works with both.
export { BoundingBox3, Frame3, Transform3, Vector3 } from 'mc-curves';

export type {
  ClosestSurfacePoint,
  CurveSurfaceIntersection,
  FirstFundamentalForm,
  ISurface,
  PrincipalCurvatures,
  SecondFundamentalForm,
  SurfaceDomain,
  SurfaceFrame,
  SurfaceIntersection,
  SurfaceIntersectionOptions,
  SurfaceJSON,
  SurfaceMesh,
  SurfaceParameter,
  SurfaceSection,
  SurfaceTessellationOptions,
} from './surfaces/ISurface.js';

export { IsoCurve, OffsetSurface, ReversedSurface, SurfaceBase, TransformedSurface, TrimmedSurface } from './surfaces/SurfaceBase.js';
export { FrameSurface } from './surfaces/FrameSurface.js';
export type { SurfaceRange } from './surfaces/FrameSurface.js';

export { ConeSurface, CylinderSurface, PlaneSurface, SphereSurface, TorusSurface } from './surfaces/Elementary.js';
export {
  EllipsoidSurface,
  EllipticParaboloidSurface,
  HelicoidSurface,
  HyperbolicParaboloidSurface,
  HyperboloidSurface,
} from './surfaces/Quadrics.js';
export { BezierSurface, BSplineSurface, NurbsSurface, revolveNurbs } from './surfaces/NurbsSurface.js';
export { BilinearSurface, GridSurface, HermiteSurface } from './surfaces/Patches.js';
export type { HermiteSurfaceData } from './surfaces/Patches.js';
export { CoonsSurface, ExtrusionSurface, GordonSurface, LoftSurface, RevolutionSurface, RuledSurface } from './surfaces/CurveSurfaces.js';
export type { GordonOptions, LoftOptions } from './surfaces/CurveSurfaces.js';
export { PipeSurface, RotationMinimizingFrames, SweepSurface } from './surfaces/SweepSurfaces.js';
export type { SweepOptions } from './surfaces/SweepSurfaces.js';
export { CompositeSurface, ParametricSurface } from './surfaces/Compound.js';
export type { ParametricSurfaceOptions } from './surfaces/Compound.js';

export { curveFromJSON, surfaceFromJSON } from './surfaces/serialization.js';
