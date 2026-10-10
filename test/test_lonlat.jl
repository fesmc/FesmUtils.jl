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
    @test Ft[1, 1] ≈ sum(F[1:4, 1:4] .* A[1:4, 1:4]) / sum(A[1:4, 1:4]) rtol=1e-5   # polar cell: sphere vs ellipsoid
    # Same source on 0:360 longitudes
    src360 = LonLatGrid(0.25:0.5:359.75, src.lat)
    F360 = circshift(F, -360)
    @test remap(tgt, src360, F360)[1] ≈ Ft
    # Target cells across the wrap of the source take cells on both sides
    shifted = LonLatGrid(-180:2:178, tgt.lat)            # first cell -181..-179
    Fs, fs = remap(shifted, src, F)
    k = [719, 720, 1, 2]
    @test Fs[1, 1] ≈ sum(F[k, 1:4] .* A[k, 1:4]) / sum(A[k, 1:4]) rtol=1e-5
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

@testset "grid angle and vector rotation" begin
    proj = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    g = ProjGrid("T", (-800, 800), (-800, 800), 100, proj)
    α = grid_angle(g)
    lon, _ = lonlat(g)
    # East is counterclockwise from the x axis by lon - lon_0 (north polar stereographic)
    @test all(abs.(mod.(α .- deg2rad.(lon .+ 45) .+ π, 2π) .- π) .< 1e-6)
    # Rotation round trip, and a northward vector points to the pole
    ue, vn = rand(size(g)...), rand(size(g)...)
    ux, uy = rotate_to_grid(α, ue, vn)
    e, n = rotate_to_geographic(α, ux, uy)
    @test e ≈ ue && n ≈ vn
    ux, uy = rotate_to_grid(α, zeros(size(g)), ones(size(g)))
    xc = [x for x in g.xc, _ in g.yc]
    yc = [y for _ in g.xc, y in g.yc]
    r = hypot.(xc, yc)
    @test all(((ux .* xc .+ uy .* yc) ./ r)[r .> 1] .≈ -1)
end

@testset "class fractions between lon-lat grids" begin
    src = LonLatGrid("S", 1.0)
    tgt = LonLatGrid("T", 2.0)
    M = Array{Union{Missing,Int}}([lo < 0 ? 1 : 2 for lo in src.lon, _ in src.lat])
    M[:, 1:10] .= missing
    fr, fv = remap_fractions(tgt, src, M, (1, 2))
    @test all(fr[1][tgt.lon .< 0, 6:end] .== 1) && all(fr[2][tgt.lon .> 0, 6:end] .== 1)
    @test all(fv[:, 1:5] .== 0) && all(fv[:, 6:end] .== 1)
end

@testset "uneven latitudes (Gaussian grid)" begin
    # Gaussian latitudes of T31 (48 roots of the Legendre polynomial), 3.75° longitudes
    function legendre_roots(n)
        map(1:n) do k
            x = cos(π * (k - 0.25) / (n + 0.5))
            for _ in 1:20
                p0, p1 = 1.0, x
                for m in 2:n
                    p0, p1 = p1, ((2m - 1) * x * p1 - (m - 1) * p0) / m
                end
                x -= p1 / (n * (x * p1 - p0) / (x^2 - 1))
            end
            x
        end
    end
    n = 48
    lat = sort(asind.(legendre_roots(n)))
    src = LonLatGrid(-180:3.75:176.25, lat)
    @test !uniform_lat(src) && src.isglobal
    @test src.latedges[1] == -90 && src.latedges[end] == 90
    @test sum(cell_area(src)) ≈ sum(cell_area(LonLatGrid("G", 5))) rtol=1e-9
    @test_throws ArgumentError LonLatGrid(0:1:9, [0.0, 2.0, 1.0])

    # Exact conservative remapping keeps the area-weighted global mean
    F = [sind(la) + cosd(lo) for lo in src.lon, la in src.lat]
    tgt = LonLatGrid("T", 2.0)
    Ft, fv = remap(tgt, src, F)
    @test all(fv .≈ 1)
    @test sum(Ft .* cell_area(tgt)) / sum(cell_area(tgt)) ≈ sum(F .* cell_area(src)) / sum(cell_area(src)) atol=1e-12

    # Sampled remapping onto a projected grid and bilinear interpolation of latitude
    proj = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    pg = ProjGrid("P", (-800, 800), (-800, 800), 50, proj)
    L = repeat(permutedims(src.lat), length(src.lon))
    Fp, fp = remap(pg, src, L; nsub=8)
    _, latc = lonlat(pg)
    @test all(fp .== 1) && maximum(abs.(Fp .- latc)) < 3      # cell means, 3.7° cells
    Fb, _ = remap_bilinear(pg, src, L)
    @test maximum(abs.(Fb .- latc)[latc .< src.lat[end]]) < 1e-6

    # Explicit edges, and a regional grid
    reg = LonLatGrid(0.5:1:9.5, [40.5, 41.4, 42.5]; latedges=[40.0, 41.0, 42.0, 43.0])
    @test reg.latedges == [40.0, 41.0, 42.0, 43.0]
    @test lat_bounds(reg) == (40.0, 43.0)
end
