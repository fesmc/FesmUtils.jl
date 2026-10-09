# FesmUtils.jl

Utilities shared across FESM components and data processing.

## Grids and remapping (implemented)

- `ProjGrid`: regular grid on a map projection (km), read from/written to cdo grid description files (`ProjGrid("grid_GRL-8KM.txt")`, `write_griddes`, `write_grid_nc`). `coarsen` and `crop` derive nested grids.
- `LonLatGrid`: regular lon-lat grid (e.g. GEBCO).
- `remap(tgt, src, F) -> (Ft, f_valid)`:
  - `src::ProjGrid` on the same projection: exact conservative remapping (separable 1D overlaps). Use `AlignedMap(tgt, src)` to reuse weights.
  - `src::LonLatGrid`: area-weighted mean from `nsub` x `nsub` samples per target cell (`nsub` keyword).
- `remap_fractions(tgt, src, M, classes)`: area fraction of each class of a categorical field.
- `distance_to(mask, dx, dy)`: exact Euclidean distance to the nearest `true` cell.

All of these are threaded: start Julia with `-t N` (or `JULIA_NUM_THREADS`).

## Planned

## FesmPreprocessing.jl

- Steps to download datasets
- Basic preprocessing to have dataset workable for our methods (data format questions, possible "general" regridding steps, but not refined to a specific domain/grid)

## FesmRemapping.jl

- functionality to take preprocessed datasets (points, polygons, grids) and map them onto the grids being used by our components (e.g. points_to_grid...). Also consider how to generate weights to be used online instead.
- Probably should use ConservativeInterpolations, Interpolations, Proj.jl, Rasters(?),....??

## FesmOperations.jl

- Any general operators that could be used in Pagos, etc.
- Sometimes use external functionality but wrap it in a convenient way. Sometimes, make a performant version (e.g. convolutions, FFTs, etc) that is specific for our fast modeling target.
- If an external package does work well, then preferred to import this directly into e.g. Pagos.jl.
