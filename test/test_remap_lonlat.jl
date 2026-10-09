@testset "remap lonlat" begin
    proj = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    tgt = ProjGrid("T", (-720, 960), (-3450, -570), 16, proj)

    # Global 0.05 deg source, northern band only
    lon = collect(-180 + 0.025 : 0.05 : 180)
    lat = collect(50 + 0.025 : 0.05 : 90)
    src = LonLatGrid(lon, lat)
    @test src.isglobal

    # Constant field is reproduced everywhere
    Ft, fv = remap(tgt, src, fill(3.0f0, size(src)...); nsub=4)
    @test all(Ft .== 3) && all(fv .== 1)

    # Latitude field: cell mean close to centre latitude
    L = repeat(Float32.(lat)', length(lon))
    Ft, fv = remap(tgt, src, L; nsub=8)
    _, latc = lonlat(tgt)
    @test maximum(abs.(Ft .- latc)) < 0.05

    # Integer input, missing values and longitude wrap
    I = fill(Int16(5), size(src)...)
    Ft, _ = remap(tgt, src, I; nsub=2)
    @test eltype(Ft) == Float32 && all(Ft .== 5)
    Fm = Array{Union{Missing,Float32}}(fill(1.0f0, size(src)...))
    Fm[:, lat .> 80] .= missing
    Ft, fv = remap(tgt, src, Fm; nsub=4)
    @test all(isnan.(Ft[latc .> 81])) && all(fv[latc .< 79] .== 1)

    # Source band not covering the domain gives zero coverage
    srcs = LonLatGrid(lon, collect(80.025:0.05:90))
    Ft, fv = remap(tgt, srcs, ones(Float32, size(srcs)...); nsub=2)
    @test all(fv[latc .< 79] .== 0)

    # Class fractions: west/east halves of the globe
    M = [l < -45 ? 1 : 2 for l in lon, _ in lat]
    fr, fv = remap_fractions(tgt, src, M, (1, 2); nsub=4)
    @test all(fr[1] .+ fr[2] .≈ 1)
    lonc, _ = lonlat(tgt)
    @test all(fr[1][lonc .< -60 .&& lonc .> -170] .== 1)

    # Same answer with any number of threads: deterministic per-cell accumulation
    @test remap(tgt, src, L; nsub=8)[1] == remap(tgt, src, L; nsub=8)[1]
end
