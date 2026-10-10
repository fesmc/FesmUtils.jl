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
