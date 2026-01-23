# FesmUtils.jl

What do we want to house here?

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
