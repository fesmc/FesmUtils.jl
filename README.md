# FesmUtils.jl

Utilities shared across FESM components and data processing.

## Grids and remapping (implemented)

- `ProjGrid`: regular grid on a map projection (km), read from/written to cdo grid description files (`ProjGrid("grid_GRL-8KM.txt")`, `write_griddes`, `write_grid_nc`). `coarsen` and `crop` derive nested grids.
- `LonLatGrid`: regular lon-lat grid (e.g. GEBCO).
- `remap(tgt, src, F) -> (Ft, f_valid)`:
  - `src::ProjGrid` on the same projection: exact conservative remapping (separable 1D overlaps). Use `AlignedMap(tgt, src)` to reuse weights.
  - `src::LonLatGrid`, or `src::ProjGrid` on any projection with the `nsub` keyword: area-weighted mean from `nsub` x `nsub` samples per target cell.
- `remap_fractions(tgt, src, M, classes)`: area fraction of each class of a categorical field (same methods).
- `xy_bounds(g, proj)`: extent of a grid in another projection (e.g. to crop the target to a source tile).
- `distance_to(mask, dx, dy)`: exact Euclidean distance to the nearest `true` cell.
- `remap_dominant(tgt, src, M)`: class of a categorical field covering the largest area of each target cell (same projection, exact overlaps).

## Masks and regions

- `rasterize(g, rings)`: cells whose centre lies inside a polygon (rings of `(x, y)` in grid coordinates, even-odd rule).
- `erode(mask, r, dx)`, `dilate(mask, r, dx)`: erosion and dilation by a distance.
- `flood_fill(mask, seeds)`: cells of `mask` connected to the seeds.
- `extend_labels(L, allowed, dx)`: nearest label along paths through `allowed` cells.
- Region codes: a region with path `"1.3.1"` (subregion 1 of region 3 of region 1) has code `region_code("1.3.1") == 10301`, two digits per level below the first. `region_level`, `region_path`, `region_ancestor`, `in_region(codes, code)` and `region_flag_attrib` (CF flag attributes) work with fields of codes.

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
