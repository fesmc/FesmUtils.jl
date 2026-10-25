# FesmUtils.jl

Utilities shared across FESM components and data processing.

## Grids and remapping (implemented)

- `ProjGrid`: regular grid on a map projection (km), read from/written to cdo grid description files (`ProjGrid("grid_GRL-8KM.txt")`, `write_griddes`, `write_grid_nc`). `coarsen` and `crop` derive nested grids. Polar stereographic (`polar_stereographic_proj`), transverse Mercator (`transverse_mercator_proj`, e.g. UTM) and Albers equal-area (`albers_equal_area_proj`, e.g. High Mountain Asia) projections; transverse Mercator and Albers grid descriptions give the projection as `proj_params` (PROJ string in m), as fesm-utils writes and reads them.
- `LonLatGrid`: regular lon-lat grid (e.g. GEBCO).
- `remap(tgt, src, F) -> (Ft, f_valid)`:
  - `src::ProjGrid` on the same projection: exact conservative remapping (separable 1D overlaps). Use `AlignedMap(tgt, src)` to reuse weights.
  - `src::LonLatGrid`, or `src::ProjGrid` on any projection with the `nsub` keyword: area-weighted mean from `nsub` x `nsub` samples per target cell.
- `remap_fractions(tgt, src, M, classes)`: area fraction of each class of a categorical field (same methods).
- `xy_bounds(g, proj)`: extent of a grid in another projection (e.g. to crop the target to a source tile).
- `distance_to(mask, dx, dy)`: exact Euclidean distance to the nearest `true` cell.
- `remap_dominant(tgt, src, M; parent)`: class of a categorical field covering the largest area of each target cell (same projection, exact overlaps), optionally nested in a coarser classification.

## Masks and regions

- `rasterize(g, rings)`, `rasterize!(mask, g, rings)`: cells whose centre lies inside a polygon (rings of `(x, y)` in grid coordinates, even-odd rule); the in-place version adds a polygon to a mask.
- `erode(mask, r, dx)`, `dilate(mask, r, dx)`: erosion and dilation by a distance.
- `flood_fill(mask, seeds)`: cells of `mask` connected to the seeds.
- `extend_labels(L, allowed, dx)`: nearest label along paths through `allowed` cells.
- `laplace_fill(F, gaps)`, `laplace_fill!(F, gaps)`: fill the gaps of a field by harmonic (Laplace) interpolation from the valid cells around them; gaps not connected to valid cells stay NaN.
- Region codes: a region with path `"1.3.1"` (subregion 1 of region 3 of region 1) has code `region_code("1.3.1") == 10301`, two digits per level below the first. `region_level`, `region_path`, `region_ancestor`, `in_region(codes, code)` and `region_flag_attrib` (CF flag attributes) work with fields of codes.

All of these are threaded: start Julia with `-t N` (or `JULIA_NUM_THREADS`).

## Namelist files (`FesmUtils.Namelist`)

Par files as fesm-utils `nml.f90` reads them (phase 0 rule: vector elements split at blanks outside quotes, so `units = "kg m-2 yr-1" "mm d-1"` is two strings; a comment after `&group` is allowed), with every comment kept:

```julia
using FesmUtils.Namelist
groups = read_namelist("yelmox/yelmox_Greenland.nml")    # Vector{Group}, Group.keys::Vector{Key}
parse_value(groups[1].keys[1].raw)                       # String, Bool, Int, Float64 or a Vector
write_namelist("copy.nml", groups)                       # read_namelist("copy.nml") == groups
```

`resolve` gives the configuration a run sees: the base par file, the groups it names (`GroupRef` rules; `FORCING_REFS` for coefficients, terms and snapshots) from the base or the library files, each group filled from its defaults group, and runme `-p` overrides, with the source of every value (`base`, `library`, `default`, `override`):

```julia
cfg = load_config(yelmox_root)                           # .runme/info.json: par_defaults, par_library
groups, order = resolve(cfg, "yelmox/yelmox_Greenland.nml", parse_overrides(["ctrl.dtt=5", "forcing.coefs==alpha,dTo"]))
value_of(groups["coef_alpha"], "method"), groups["coef_alpha"].entries["method"].source
```

`Config(root; schema, library)` sets the defaults groups (`name or Regex => (file, group)`) and the library explicitly. The module docstring has the parsing rules and what is not supported (flat `group.key = value` lines).

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
