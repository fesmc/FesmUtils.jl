@testset "remap conservative" begin
    utm = transverse_mercator_proj(lat_0=0, lon_0=-75, k=0.9996, x_0=500000, y_0=10000000,
                                   a=6378137, rf=298.257223563)
    # 2 km of 25 m cells near 54.4S, 70.8W, and a 1 arc second source over it
    tgt = ProjGrid("T", (766.0125, 767.9875), (3962.0125, 3963.9875), 0.025, utm)
    d = 1 / 3600
    lon = collect(-70.95 + d / 2 : d : -70.8)
    lat = collect(-54.45 + d / 2 : d : -54.35)
    src = LonLatGrid(lon, lat)
    m = ConservativeMap(tgt, src)

    # Weights cover each target cell exactly
    Ft, fv = remap(m, ones(Float32, size(src)...))
    @test all(Ft .== 1) && all(isapprox.(fv, 1; atol=1e-9))

    # A linear field in the projected plane has the cell-centre value as its cell mean
    trans = FesmUtils.Proj.Transformation("EPSG:4326", utm; always_xy=true)
    xy = [trans((lo, la)) for lo in lon, la in lat]
    Fx = first.(xy)
    Ft, _ = remap(m, Fx)
    xc = [x for x in tgt.xc, _ in tgt.yc]
    @test maximum(abs.(Ft .- xc)) < 2e-3   # km; source cells of about 30 m

    # Agrees with dense sampling, and has no aliasing: the spread of a field that varies
    # at the source scale is smooth across cells
    Fr = Float64[sin(2pi * i / 7) + cos(2pi * j / 5) for i in eachindex(lon), j in eachindex(lat)]
    Fc, _ = remap(m, Fr)
    Fs, _ = remap(tgt, src, Fr; nsub=64)
    @test maximum(abs.(Fc .- Fs)) < 0.02

    # Missing values are excluded, with the valid area fraction
    Fm = Array{Union{Missing,Float32}}(fill(2.0f0, size(src)...))
    Fm[:, lat .> -54.41] .= missing
    Ft, fv = remap(m, Fm)
    @test all(v -> isnan(v) || v == 2, Ft)
    @test all(f -> -1e-9 <= f <= 1 + 1e-9, fv) && any(f -> 0 < f < 1, fv)

    # Target cells outside the source have no valid area
    srcs = LonLatGrid(lon, collect(-54.45 + d / 2 : d : -54.41))
    Ft, fv = remap(ConservativeMap(tgt, srcs), ones(Float32, size(srcs)...))
    _, latc = lonlat(tgt)
    @test count(latc .> -54.409) > 0
    @test all(fv[latc .> -54.409] .== 0) && all(isnan.(Ft[latc .> -54.409]))

    # Source cells larger than the target cells
    src4 = LonLatGrid(collect(-71 + 0.05 : 0.1 : -70.6), collect(-54.6 + 0.05 : 0.1 : -54.2))
    F4 = Float64[i + 10j for i in 1:4, j in 1:4]
    Ft, fv = remap(ConservativeMap(tgt, src4), F4)
    Fs, _ = remap(tgt, src4, F4; nsub=64)
    @test all(isapprox.(fv, 1; atol=1e-9)) && maximum(abs.(Ft .- Fs)) < 0.2
    @test_throws DimensionMismatch remap(m, ones(3, 3))
end
