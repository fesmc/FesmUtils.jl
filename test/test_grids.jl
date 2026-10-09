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
end
