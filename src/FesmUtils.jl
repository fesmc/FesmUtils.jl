module FesmUtils

using NCDatasets
import Proj

export ProjGrid, LonLatGrid, spacing, grid_name, coarsen, crop
export polar_stereographic_proj, cf_grid_mapping
export write_griddes, write_grid_nc, init_grid_nc!, lonlat, cell_area, lat_bounds

include("grids.jl")

end # module
