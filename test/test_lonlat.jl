@testset "lon-lat grids" begin
    g = LonLatGrid("GLOBAL-0.5DEG", 0.5)
    @test size(g) == (720, 360) && g.isglobal
    @test g.lon[1] == -179.75 && g.lat[end] == 89.75
    @test grid_name("GLOBAL", 0.25; units="deg") == "GLOBAL-0.25DEG"
    @test grid_name("GLOBAL", 2.0; units="deg") == "GLOBAL-2DEG"
    @test_throws ArgumentError LonLatGrid("X", 0.7)
    @test grid_dims(g) == ("lon", "lat")
    @test lat_bounds(g) == (-90.0, 90.0)

    # Cell areas sum to the area of the WGS84 ellipsoid
    @test sum(cell_area(LonLatGrid("G", 5))) ≈ 5.10065621724e14 rtol=1e-9

    # Grid files
    mktempdir() do dir
        write_griddes(joinpath(dir, "g.txt"), g)
        @test occursin("gridtype = lonlat", read(joinpath(dir, "g.txt"), String))
        write_grid_nc(joinpath(dir, "g.nc"), g)
        NCDataset(joinpath(dir, "g.nc")) do ds
            @test ds.attrib["grid_name"] == "GLOBAL-0.5DEG"
            @test size(ds["area"]) == (720, 360)
        end
    end

    # Conservative lon-lat to lon-lat: exact, weighted by area on the sphere (within
    # 1e-7 of the ellipsoid), any longitude range
    src = LonLatGrid("S", 0.5)
    F = rand(size(src)...)
    tgt = LonLatGrid("T", 2)
    Ft, fv = remap(tgt, src, F)
    @test all(fv .≈ 1)
    A, At = cell_area(src), cell_area(tgt)
    @test sum(Ft .* At) ≈ sum(F .* A) rtol=1e-6
    @test Ft[1, 1] ≈ sum(F[1:4, 1:4] .* A[1:4, 1:4]) / sum(A[1:4, 1:4]) rtol=1e-6
    # Same source on 0:360 longitudes
    src360 = LonLatGrid(0.25:0.5:359.75, src.lat)
    F360 = circshift(F, -360)
    @test remap(tgt, src360, F360)[1] ≈ Ft
    # Target cells across the wrap of the source take cells on both sides
    shifted = LonLatGrid(-180:2:178, tgt.lat)            # first cell -181..-179
    Fs, fs = remap(shifted, src, F)
    k = [719, 720, 1, 2]
    @test Fs[1, 1] ≈ sum(F[k, 1:4] .* A[k, 1:4]) / sum(A[k, 1:4]) rtol=1e-6
    @test all(fs .≈ 1)
    # Regional target and source: only the overlap counts
    reg = LonLatGrid(10.25:0.5:19.75, 40.25:0.5:49.75)
    Fr, fr = remap(LonLatGrid(5:2:23, 37:2:53), reg, ones(size(reg)))
    @test all(Fr[fr .> 0] .≈ 1) && fr[1, 1] == 0 && fr[5, 5] ≈ 1

    # Projected source onto a lon-lat target (sampled)
    proj = polar_stereographic_proj(lat_0=-90, lat_ts=-71, lon_0=0, a=6378137, rf=298.257223563)
    pg = ProjGrid("P", (-3000, 3000), (-3000, 3000), 20, proj)
    lon, lat = lonlat(pg)
    Fl, fl = remap(LonLatGrid("G", 1), pg, lat; nsub=8)
    ok = fl .> 0.999
    @test all(abs.(Fl[ok] .- repeat(permutedims(LonLatGrid("G", 1).lat), 360)[ok]) .< 0.6)
    @test all(fl[:, 100:end] .== 0)
end
