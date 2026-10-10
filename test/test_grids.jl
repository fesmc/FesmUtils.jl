const PS_NORTH = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)

@testset "grids" begin
    g = ProjGrid("GRL-8KM", (-720, 960), (-3450, -570), 8, PS_NORTH)
    @test size(g) == (211, 361)
    @test spacing(g) == (8.0, 8.0)

    @test_throws ArgumentError ProjGrid("bad", [0.0, 1.0, 3.0], [0.0, 1.0], PS_NORTH)
    @test_throws ArgumentError ProjGrid("bad", [0.0, 1.0], [0.0, 1.0], replace(PS_NORTH, "km" => "m"))

    @test grid_name("GRL", 0.5) == "GRL-500M"
    @test grid_name("ANT", 1) == "ANT-1KM"
    @test grid_name("NH", 32.0) == "NH-32KM"

    # Coarse cell centres are a subset of fine cell centres
    c = coarsen(g, 4; name="GRL-32KM")
    @test size(c) == (53, 91)
    @test all(x -> x in g.xc, c.xc) && all(y -> y in g.yc, c.yc)

    s = crop(g, (-400, 400), (-3050, -1050); name="sub")
    @test s.xc[1] == -400 && s.xc[end] == 400 && s.yc[1] == -3050 && s.yc[end] == -1050

    # griddes round trip
    path = tempname() * ".txt"
    write_griddes(path, g)
    g2 = ProjGrid(path; name="GRL-8KM")
    @test g2.xc == g.xc && g2.yc == g.yc
    @test cf_grid_mapping(g2) == cf_grid_mapping(g)

    # Areal scale factor is 1 at the standard parallel, so cells there are dx*dy
    lon, lat = lonlat(g)
    area = cell_area(g, lon, lat)
    k = argmin(abs.(lat .- 70))
    @test isapprox(area[k], 64e6; rtol=1e-3)

    latmin, latmax = lat_bounds(g; margin=0)
    @test latmin <= minimum(lat) && latmax >= maximum(lat)
    gs = ProjGrid("ANT", (-3040, 3040), (-3040, 3040), 32,
                  polar_stereographic_proj(lat_0=-90, lat_ts=-71, lon_0=0, a=6378137, rf=298.257223563))
    @test lat_bounds(gs)[1] == -90.0

    ncpath = tempname() * ".nc"
    write_grid_nc(ncpath, g)
    @test isfile(ncpath)
    NCDataset(ncpath) do ds
        @test eltype(ds["area"].var) == Float32
        @test maximum(abs.(ds["area"][:, :] .- area) ./ area) < 1e-6
        @test maximum(abs.(ds["lat2D"][:, :] .- lat)) < 1e-4
    end

    # Transverse Mercator (UTM 18S, the SRG grids of yelmox)
    UTM18S = transverse_mercator_proj(lat_0=0, lon_0=-75, k=0.9996, x_0=500000, y_0=10000000,
                                      a=6378137, rf=298.257223563)
    t = ProjGrid("SRG-16KM", (576.962171469987, 640.962171469987), (4809.977314639904, 4841.977314639904),
                 16, UTM18S)
    @test size(t) == (5, 3)
    cf = Dict(cf_grid_mapping(t))
    @test cf["grid_mapping_name"] == "transverse_mercator"
    @test cf["longitude_of_central_meridian"] == -75 && cf["latitude_of_projection_origin"] == 0
    @test cf["scale_factor_at_central_meridian"] == 0.9996
    @test cf["false_easting"] == 500 && cf["false_northing"] == 10000

    # Cell centres in km match EPSG:32718 in m
    lon, lat = lonlat(t)
    utm = FesmUtils.Proj.Transformation("EPSG:32718", "EPSG:4326"; always_xy=true)
    lo, la = utm((t.xc[2] * 1000, t.yc[3] * 1000))
    @test isapprox(lon[2, 3], lo; atol=1e-9) && isapprox(lat[2, 3], la; atol=1e-9)

    # griddes: proj_params in m, as fesm-utils writes and reads it, and round trip
    path = tempname() * ".txt"
    write_griddes(path, t)
    txt = read(path, String)
    @test occursin("proj_params = \"+proj=tmerc +lon_0=-75 +lat_0=0 +k=0.9996 +x_0=500000 " *
                   "+y_0=10000000 +a=6378137 +rf=298.257223563 +units=m\"", txt)
    @test !occursin("grid_mapping_name", txt)
    t2 = ProjGrid(path; name="SRG-16KM")
    @test t2.xc ≈ t.xc && t2.yc ≈ t.yc && same_projection(t2, t)
    @test cf_grid_mapping(t2) == cf_grid_mapping(t)

    # A grid description written by fesm-utils (v1 maps/grid_SRG-16KM.txt)
    path = tempname() * ".txt"
    write(path, """
        # fesmutils grid description
        # mtype  = transverse_mercator
        # units  = kilometers
        gridtype = projection
        gridsize =         15
        xsize    =          5
        ysize    =          3
        xname    = xc
        xunits   = km
        yname    = yc
        yunits   = km
        xfirst   = 576.962171469987
        xinc     = 16.000000000000
        yfirst   = 4809.977314639904
        yinc     = 16.000000000000
        grid_mapping = crs
        proj_params = "+proj=tmerc +lon_0=-75.000000000 +lat_0=.000000000 +k=.999600000 +x_0=500000.000000 +y_0=10000000.000000 +a=6378137.000 +rf=298.257223563 +units=m"
        """)
    t3 = ProjGrid(path; name="SRG-16KM")
    @test size(t3) == (5, 3) && t3.xc ≈ t.xc && t3.yc ≈ t.yc && same_projection(t3, t)

    ncpath = tempname() * ".nc"
    write_grid_nc(ncpath, t)
    NCDataset(ncpath) do ds
        @test ds["crs"].attrib["grid_mapping_name"] == "transverse_mercator"
        @test ds["crs"].attrib["proj_params"] == UTM18S
        @test maximum(abs.(ds["lat2D"][:, :] .- lat)) < 1e-4
    end
end
