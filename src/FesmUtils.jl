module FesmUtils

using NCDatasets
import Proj

export ProjGrid, LonLatGrid, spacing, grid_name, coarsen, crop
export polar_stereographic_proj, cf_grid_mapping
export write_griddes, write_grid_nc, init_grid_nc!, lonlat, cell_area, lat_bounds
export AlignedMap, same_projection, remap, remap_fractions

include("grids.jl")
include("remap_aligned.jl")
include("remap_lonlat.jl")

end # module
