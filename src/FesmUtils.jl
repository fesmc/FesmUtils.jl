module FesmUtils

using NCDatasets
import Proj
using SparseArrays: sparse

export ProjGrid, LonLatGrid, spacing, grid_name, coarsen, crop
export polar_stereographic_proj, transverse_mercator_proj, albers_equal_area_proj, cf_grid_mapping
export write_griddes, write_grid_nc, init_grid_nc!, grid_dims, lonlat, cell_area, lat_bounds, xy_bounds
export grid_angle, rotate_to_grid, rotate_to_geographic, uniform_lat, gaussian_grid
export AlignedMap, SampledMap, same_projection, remap, remap_fractions, remap_dominant, remap_bilinear
export smooth
export distance_to, erode, dilate, flood_fill, extend_labels
export laplace_fill, laplace_fill!
export rasterize, rasterize!
export region_code, region_level, region_path, region_ancestor, in_region, region_flag_attrib

include("grids.jl")
include("remap_aligned.jl")
include("remap_sampled.jl")
include("remap_bilinear.jl")
include("smooth.jl")
include("distance.jl")
include("morphology.jl")
include("fill.jl")
include("rasterize.jl")
include("region_codes.jl")

end # module
