@testset "remap bilinear" begin
    # Lon-lat to lon-lat: linear fields are reproduced, periodic in longitude
    src = LonLatGrid("S", 1)
    tgt = LonLatGrid("T", 0.25)
    lat = repeat(permutedims(src.lat), length(src.lon))
    Ft, fv = remap_bilinear(tgt, src, lat)
    tl = repeat(permutedims(tgt.lat), length(tgt.lon))
    inner = abs.(tl) .< 89.5
    @test all(Ft[inner] .≈ tl[inner]) && all(fv .≈ 1)
    @test all(Ft[.!inner] .≈ clamp.(tl[.!inner], -89.5, 89.5))     # nearest row at the poles
    L = zeros(size(src)); L[1, :] .= 1                               # cell at -179.5
    Fl, _ = remap_bilinear(tgt, src, L)
    @test Fl[end, 90] ≈ 0.375 rtol=1e-12             # 179.875: 3/8 of the way from 179.5 to 180.5
    @test Fl[1, 90] ≈ 0.625 rtol=1e-12               # -179.875

    # Projected to projected (same projection): linear in x and y
    proj = polar_stereographic_proj(lat_0=-90, lat_ts=-71, lon_0=0, a=6378137, rf=298.257223563)
    ps = ProjGrid("PS", (-1000, 1000), (-1000, 1000), 20, proj)
    pt = ProjGrid("PT", (-993, 993), (-993, 993), 6, proj)
    X = [x + 2y for x in ps.xc, y in ps.yc]
    Xt, xv = remap_bilinear(pt, ps, X)
    @test all(Xt .≈ [x + 2y for x in pt.xc, y in pt.yc]) && all(xv .≈ 1)

    # Missing neighbours are left out
    Xm = copy(X); Xm[51, 51] = NaN
    Xmt, mv = remap_bilinear(pt, ps, Xm)
    @test all(isfinite, Xmt) && minimum(mv) < 1

    # Lon-lat to projected, and outside the source
    lat2 = lonlat(pt)[2]
    Lt, lv = remap_bilinear(pt, src, lat)
    @test maximum(abs, (Lt .- lat2)[lat2 .> -89]) < 0.01      # nearest row within 0.5° of the pole
    reg = LonLatGrid(0.5:1:9.5, -49.5:1:-40.5)
    _, rv = remap_bilinear(pt, reg, ones(size(reg)))
    @test all(rv .== 0)
end

@testset "remap sd" begin
    # Projected source, target cells smaller than the source cells: a linear field has
    # the spread of the samples of a uniform slope in every cell (no aliasing)
    proj = polar_stereographic_proj(lat_0=-90, lat_ts=-71, lon_0=0, a=6378137, rf=298.257223563)
    ps = ProjGrid("PS", (-100, 100), (-100, 100), 4, proj)
    pt = ProjGrid("PT", (-49.5, 49.5), (-49.5, 49.5), 1, proj)
    X = [3x for x in ps.xc, y in ps.yc]
    nsub = 8
    Xm, sd, fv = remap_sd(pt, ps, X; nsub)
    @test all(sd .≈ 3 * sqrt((nsub^2 - 1) / 12) / nsub) && all(fv .== 1)
    @test all(Xm .≈ [3x for x in pt.xc, y in pt.yc])
    Cm, sd, _ = remap_sd(pt, ps, fill(5.0f0, size(ps)...); nsub=4)
    @test all(sd .== 0) && all(Cm .== 5) && eltype(sd) == Float32

    # Lon-lat source: a field linear in latitude, sd close to its slope times the cell
    # size over sqrt(12)
    utm = transverse_mercator_proj(lat_0=0, lon_0=-75, k=0.9996, x_0=500000, y_0=10000000,
                                   a=6378137, rf=298.257223563)
    tgt = ProjGrid("T", (766.0125, 767.9875), (3962.0125, 3963.9875), 0.025, utm)
    d = 1 / 3600
    src = LonLatGrid(collect(-70.95 + d / 2 : d : -70.8), collect(-54.45 + d / 2 : d : -54.35))
    L = repeat(permutedims(src.lat .* 111195), length(src.lon))    # metres northward
    _, sd, fv = remap_sd(tgt, src, L; nsub=16)
    @test all(fv .== 1) && all(abs.(sd ./ (25 / sqrt(12)) .- 1) .< 0.01)

    # Outside the source, and missing values
    Lm = copy(L); Lm[:, src.lat .> -54.41] .= NaN
    _, sd, fv = remap_sd(tgt, src, Lm; nsub=4)
    _, latc = lonlat(tgt)
    @test all(fv[latc .> -54.405] .== 0) && all(isnan.(sd[latc .> -54.405])) && all(fv[latc .< -54.415] .== 1)
    @test_throws DimensionMismatch remap_sd(tgt, src, ones(3, 3); nsub=2)
end
