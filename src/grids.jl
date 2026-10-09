# Grid definitions: regular projected grids (km) and regular lon-lat grids (degrees).

const AXIS_RTOL = 1e-6

function _check_axis(v::AbstractVector, name::AbstractString)
    length(v) >= 2 || throw(ArgumentError("$name must have at least 2 points"))
    d = v[2] - v[1]
    d > 0 || throw(ArgumentError("$name must be ascending"))
    for k in 2:length(v)
        isapprox(v[k] - v[k-1], d; rtol=AXIS_RTOL) ||
            throw(ArgumentError("$name must be uniformly spaced"))
    end
    return nothing
end

"""
    ProjGrid(name, xc, yc, proj)
    ProjGrid(name, xlim, ylim, dx, proj)

Regular grid on a map projection. `xc` and `yc` are ascending, uniformly spaced
cell centres in km, and `proj` is a PROJ string with `+units=km`.
"""
struct ProjGrid
    name::String
    xc::Vector{Float64}
    yc::Vector{Float64}
    proj::String

    function ProjGrid(name::AbstractString, xc::AbstractVector, yc::AbstractVector, proj::AbstractString)
        _check_axis(xc, "xc")
        _check_axis(yc, "yc")
        _proj_params(proj)["units"] == "km" ||
            throw(ArgumentError("proj string must use +units=km: $proj"))
        return new(String(name), collect(Float64, xc), collect(Float64, yc), String(proj))
    end
end

function ProjGrid(name::AbstractString, xlim, ylim, dx::Real, proj::AbstractString)
    return ProjGrid(name, _axis(xlim, dx), _axis(ylim, dx), proj)
end

function _axis(lim, d)
    n = (lim[2] - lim[1]) / d
    isapprox(n, round(n); atol=1e-6) ||
        throw(ArgumentError("extent $lim is not a multiple of spacing $d"))
    return lim[1] .+ d .* (0:round(Int, n))
end

Base.size(g::ProjGrid) = (length(g.xc), length(g.yc))
spacing(g::ProjGrid) = (g.xc[2] - g.xc[1], g.yc[2] - g.yc[1])

"""
    LonLatGrid(lon, lat)

Regular lon-lat grid with ascending, uniformly spaced cell centres in degrees.
"""
struct LonLatGrid
    lon::Vector{Float64}
    lat::Vector{Float64}
    dlon::Float64
    dlat::Float64
    isglobal::Bool

    function LonLatGrid(lon::AbstractVector, lat::AbstractVector)
        _check_axis(lon, "lon")
        _check_axis(lat, "lat")
        dlon = (lon[end] - lon[1]) / (length(lon) - 1)
        dlat = (lat[end] - lat[1]) / (length(lat) - 1)
        isglobal = isapprox(length(lon) * dlon, 360.0; rtol=AXIS_RTOL)
        return new(collect(Float64, lon), collect(Float64, lat), dlon, dlat, isglobal)
    end
end

Base.size(g::LonLatGrid) = (length(g.lon), length(g.lat))

# ---------------------------------------------------------------------------
# Naming, coarsening and cropping
# ---------------------------------------------------------------------------

"""
    grid_name(prefix, dx)

Standard grid name for resolution `dx` (km): `grid_name("GRL", 0.5) == "GRL-500M"`,
`grid_name("ANT", 8) == "ANT-8KM"`.
"""
function grid_name(prefix::AbstractString, dx::Real)
    if dx < 1
        m = dx * 1000
        isinteger(round(m; digits=6)) || throw(ArgumentError("dx=$dx km is not a whole number of metres"))
        return "$prefix-$(round(Int, m))M"
    end
    isinteger(dx) || throw(ArgumentError("dx=$dx km is not a whole number of km"))
    return "$prefix-$(Int(dx))KM"
end

"""
    coarsen(g, factor; name)

Grid with spacing `factor` times that of `g`, sharing its first cell centre, so
that every coarse cell centre is also a cell centre of `g`.
"""
function coarsen(g::ProjGrid, factor::Integer; name::AbstractString)
    factor >= 1 || throw(ArgumentError("factor must be >= 1"))
    dx, dy = spacing(g)
    nx = (length(g.xc) - 1) ÷ factor + 1
    ny = (length(g.yc) - 1) ÷ factor + 1
    xc = g.xc[1] .+ (factor * dx) .* (0:nx-1)
    yc = g.yc[1] .+ (factor * dy) .* (0:ny-1)
    return ProjGrid(name, xc, yc, g.proj)
end

"""
    crop(g, xlim, ylim; name)

Subgrid of `g` with cell centres inside `xlim` and `ylim` (km, inclusive).
"""
function crop(g::ProjGrid, xlim, ylim; name::AbstractString)
    dx, dy = spacing(g)
    ix = findall(x -> xlim[1] - 1e-6dx <= x <= xlim[2] + 1e-6dx, g.xc)
    iy = findall(y -> ylim[1] - 1e-6dy <= y <= ylim[2] + 1e-6dy, g.yc)
    (isempty(ix) || isempty(iy)) && throw(ArgumentError("crop of $(g.name) is empty"))
    return ProjGrid(name, g.xc[ix], g.yc[iy], g.proj)
end

# ---------------------------------------------------------------------------
# PROJ strings and CF grid mappings
# ---------------------------------------------------------------------------

function _proj_params(proj::AbstractString)
    p = Dict{String,String}()
    for tok in split(proj)
        startswith(tok, "+") || continue
        kv = split(tok[2:end], "="; limit=2)
        p[kv[1]] = length(kv) == 2 ? kv[2] : ""
    end
    return p
end

_float(s::AbstractString) = parse(Float64, s)

function _num(s::AbstractString)
    v = parse(Float64, s)
    return isinteger(v) ? Int(v) : v
end

"""
    polar_stereographic_proj(; lat_0, lat_ts, lon_0, a, rf, x_0=0, y_0=0)

PROJ string (km units) for a polar stereographic projection. False easting and
northing are in m, as in PROJ.
"""
function polar_stereographic_proj(; lat_0, lat_ts, lon_0, a, rf, x_0=0, y_0=0)
    return "+proj=stere +lat_0=$lat_0 +lat_ts=$lat_ts +lon_0=$lon_0 " *
           "+x_0=$x_0 +y_0=$y_0 +a=$a +rf=$rf +units=km"
end

"""
    cf_grid_mapping(g)

CF grid-mapping attributes of `g` as a vector of pairs (false easting/northing in km).
"""
function cf_grid_mapping(g::ProjGrid)
    p = _proj_params(g.proj)
    p["proj"] == "stere" || error("only polar stereographic grids are supported, got +proj=$(p["proj"])")
    return [
        "grid_mapping_name" => "polar_stereographic",
        "straight_vertical_longitude_from_pole" => _float(p["lon_0"]),
        "latitude_of_projection_origin" => _float(p["lat_0"]),
        "standard_parallel" => _float(p["lat_ts"]),
        "false_easting" => _float(get(p, "x_0", "0")) / 1000,
        "false_northing" => _float(get(p, "y_0", "0")) / 1000,
        "semi_major_axis" => _float(p["a"]),
        "inverse_flattening" => _float(p["rf"]),
    ]
end

# ---------------------------------------------------------------------------
# cdo grid description files
# ---------------------------------------------------------------------------

"""
    ProjGrid(griddes; name)

Read a cdo grid description file (`grid_<NAME>.txt`). The name defaults to the
`<NAME>` part of the filename.
"""
function ProjGrid(griddes::AbstractString;
                  name::AbstractString=replace(basename(griddes), r"^grid_" => "", r"\.txt$" => ""))
    d = Dict{String,String}()
    for line in eachline(griddes)
        line = strip(line)
        (isempty(line) || startswith(line, "#") || !occursin("=", line)) && continue
        k, v = strip.(split(line, "="; limit=2))
        d[k] = v
    end
    d["gridtype"] == "projection" || error("$griddes: gridtype must be projection")
    d["grid_mapping_name"] == "polar_stereographic" || error("$griddes: only polar_stereographic is supported")
    (d["xunits"] == "km" && d["yunits"] == "km") || error("$griddes: axis units must be km")

    proj = polar_stereographic_proj(
        lat_0  = _num(d["latitude_of_projection_origin"]),
        lat_ts = _num(d["standard_parallel"]),
        lon_0  = _num(d["straight_vertical_longitude_from_pole"]),
        a      = _num(d["semi_major_axis"]),
        rf     = _num(d["inverse_flattening"]),
        x_0    = _num(get(d, "false_easting", "0")) * 1000,
        y_0    = _num(get(d, "false_northing", "0")) * 1000,
    )
    nx, ny = parse(Int, d["xsize"]), parse(Int, d["ysize"])
    x0, dx = parse(Float64, d["xfirst"]), parse(Float64, d["xinc"])
    y0, dy = parse(Float64, d["yfirst"]), parse(Float64, d["yinc"])
    return ProjGrid(name, x0 .+ dx .* (0:nx-1), y0 .+ dy .* (0:ny-1), proj)
end

"""
    write_griddes(path, g)

Write the cdo grid description file of `g`.
"""
function write_griddes(path::AbstractString, g::ProjGrid)
    nx, ny = size(g)
    dx, dy = spacing(g)
    f(v) = v isa Integer ? string(v) : string(round(v; digits=8))
    open(path, "w") do io
        println(io, "gridtype = projection")
        println(io, "gridsize = $(nx*ny)")
        println(io, "xsize    = $nx")
        println(io, "ysize    = $ny")
        println(io, "xname    = xc")
        println(io, "xunits   = km")
        println(io, "yname    = yc")
        println(io, "yunits   = km")
        println(io, "xfirst   = $(f(g.xc[1]))")
        println(io, "xinc     = $(f(dx))")
        println(io, "yfirst   = $(f(g.yc[1]))")
        println(io, "yinc     = $(f(dy))")
        println(io, "grid_mapping = crs")
        for (k, v) in cf_grid_mapping(g)
            println(io, "$k = $v")
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# Geographic coordinates and cell areas
# ---------------------------------------------------------------------------

# Per-task inverse projection (Proj transformations are not thread-safe).
_inverse_transform(g::ProjGrid) =
    Proj.Transformation(g.proj, "EPSG:4326"; always_xy=true, ctx=Proj.proj_context_create())

# Split 1:n into chunks for threaded loops.
_chunks(n::Integer) = Iterators.partition(1:n, max(1, cld(n, 4 * Threads.nthreads())))

"""
    lonlat(g) -> (lon2D, lat2D)

Longitude and latitude of the cell centres of `g` (threaded).
"""
function lonlat(g::ProjGrid)
    nx, ny = size(g)
    lon = Matrix{Float64}(undef, nx, ny)
    lat = Matrix{Float64}(undef, nx, ny)
    tasks = map(_chunks(ny)) do js
        Threads.@spawn begin
            trans = _inverse_transform(g)
            for j in js, i in 1:nx
                lon[i, j], lat[i, j] = trans((g.xc[i], g.yc[j]))
            end
        end
    end
    foreach(wait, tasks)
    return lon, lat
end

"""
    cell_area(g, lon2D, lat2D)

Ellipsoidal area (m²) of each cell of `g`: projected area divided by the areal
scale factor of the projection at the cell centre (threaded).
"""
function cell_area(g::ProjGrid, lon::AbstractMatrix, lat::AbstractMatrix)
    dx, dy = spacing(g)
    nx, ny = size(g)
    area = Matrix{Float64}(undef, nx, ny)
    tasks = map(_chunks(ny)) do js
        Threads.@spawn begin
            ctx = Proj.proj_context_create()
            P = Proj.proj_create(g.proj, ctx)
            for j in js, i in 1:nx
                fac = Proj.proj_factors(P, Proj.Coord(deg2rad(lon[i, j]), deg2rad(lat[i, j])))
                area[i, j] = dx * dy * 1e6 / fac.areal_scale
            end
        end
    end
    foreach(wait, tasks)
    return area
end

"""
    lat_bounds(g; margin=0.1)

Latitude range (degrees) covered by `g`, widened by `margin`. Used to read only the
needed latitude band of a global lon-lat dataset.
"""
function lat_bounds(g::ProjGrid; margin::Real=0.1)
    dx, dy = spacing(g)
    xe = vcat(g.xc .- dx / 2, g.xc[end] + dx / 2)
    ye = vcat(g.yc .- dy / 2, g.yc[end] + dy / 2)
    trans = _inverse_transform(g)
    lats = Float64[]
    for x in xe, y in (ye[1], ye[end])
        push!(lats, trans((x, y))[2])
    end
    for y in ye, x in (xe[1], xe[end])
        push!(lats, trans((x, y))[2])
    end
    latmin, latmax = extrema(lats)

    # A pole inside the domain is not on the boundary
    fwd = Proj.Transformation("EPSG:4326", g.proj; always_xy=true)
    for pole in (90.0, -90.0)
        xp, yp = fwd((0.0, pole))
        if xe[1] <= xp <= xe[end] && ye[1] <= yp <= ye[end]
            pole > 0 ? (latmax = 90.0) : (latmin = -90.0)
        end
    end
    return (max(latmin - margin, -90.0), min(latmax + margin, 90.0))
end

"""
    xy_bounds(g, proj; npts=16) -> ((xmin, xmax), (ymin, ymax))

Extent (km) of the cells of `g` in the projection `proj` (a PROJ string with
`+units=km`), from `npts` points along each side of its outline. Exact for conformal
projections, whose coordinates take their extremes on the outline.
"""
function xy_bounds(g::ProjGrid, proj::AbstractString; npts::Integer=16)
    dx, dy = spacing(g)
    x0, x1 = g.xc[1] - dx / 2, g.xc[end] + dx / 2
    y0, y1 = g.yc[1] - dy / 2, g.yc[end] + dy / 2
    trans = Proj.Transformation(g.proj, proj; always_xy=true)
    s = range(0, 1; length=npts)
    pts = vcat([(x0 + t * (x1 - x0), y) for t in s, y in (y0, y1)][:],
               [(x, y0 + t * (y1 - y0)) for t in s, x in (x0, x1)][:])
    xy = map(trans, pts)
    return extrema(first.(xy)), extrema(last.(xy))
end

# ---------------------------------------------------------------------------
# NetCDF grid file
# ---------------------------------------------------------------------------

"""
    write_grid_nc(path, g)

Write a NetCDF grid file with `xc`, `yc`, `crs`, `lon2D`, `lat2D` and `area`.
"""
function write_grid_nc(path::AbstractString, g::ProjGrid)
    lon, lat = lonlat(g)
    area = cell_area(g, lon, lat)
    NCDataset(path, "c") do ds
        init_grid_nc!(ds, g)
        defVar(ds, "lon2D", lon, ("xc", "yc"); attrib=["units" => "degrees_east", "grid_mapping" => "crs"])
        defVar(ds, "lat2D", lat, ("xc", "yc"); attrib=["units" => "degrees_north", "grid_mapping" => "crs"])
        defVar(ds, "area", area, ("xc", "yc"); attrib=["units" => "m^2", "grid_mapping" => "crs"])
    end
    return path
end

"""
    init_grid_nc!(ds, g)

Define the `xc`, `yc` axes and the `crs` grid-mapping variable of `g` in an open dataset.
"""
function init_grid_nc!(ds::NCDataset, g::ProjGrid)
    ds.attrib["grid_name"] = g.name
    defVar(ds, "xc", g.xc, ("xc",); attrib=[
        "standard_name" => "projection_x_coordinate", "units" => "km", "axis" => "X"])
    defVar(ds, "yc", g.yc, ("yc",); attrib=[
        "standard_name" => "projection_y_coordinate", "units" => "km", "axis" => "Y"])
    crs = defVar(ds, "crs", Int32, (); attrib=vcat(cf_grid_mapping(g), ["proj_params" => g.proj]))
    crs[] = Int32(0)
    return ds
end
