@testset "remap sampled" begin
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

@testset "remap sampled projected" begin
    ps = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    utm = "+proj=utm +zone=28 +datum=WGS84 +units=km +no_defs"

    # 1 km UTM source around the central meridian (Iceland), 4 km polar stereographic target
    src = ProjGrid("U", (400.5, 599.5), (7000.5, 7199.5), 1, utm)
    (xl, yl) = xy_bounds(src, ps)
    tgt = ProjGrid("T", (4 * floor(xl[1] / 4) - 8, 4 * ceil(xl[2] / 4) + 8),
                   (4 * floor(yl[1] / 4) - 8, 4 * ceil(yl[2] / 4) + 8), 4, ps)

    # Extent in the target projection contains all source cells
    fwd = FesmUtils.Proj.Transformation(utm, ps; always_xy=true)
    xy = [fwd((x, y)) for x in src.xc, y in src.yc]
    @test all(xl[1] .< first.(xy) .< xl[2]) && all(yl[1] .< last.(xy) .< yl[2])
    @test xy_bounds(tgt, ps) == ((tgt.xc[1] - 2, tgt.xc[end] + 2), (tgt.yc[1] - 2, tgt.yc[end] + 2))

    # Linear field: cell mean close to the value at the cell centre
    X = Float32.(repeat(src.xc, 1, length(src.yc)))
    Ft, fv = remap(tgt, src, X; nsub=8)
    inv = FesmUtils.Proj.Transformation(ps, utm; always_xy=true)
    full = fv .== 1
    xu = [inv((x, y))[1] for x in tgt.xc, y in tgt.yc]
    @test count(full) > 1000
    @test maximum(abs.(Ft[full] .- xu[full])) < 0.05
    @test all(isnan.(Ft[fv .== 0]))

    # Disk of valid data: covered area conserved (true areas)
    D = [hypot(x - 500, y - 7100) < 80 ? 1.0f0 : NaN32 for x in src.xc, y in src.yc]
    Ft, fv = remap(tgt, src, D; nsub=8)
    lon, lat = lonlat(tgt)
    area = cell_area(tgt, lon, lat) ./ 1e6
    @test sum(fv .* area) ≈ π * 80^2 rtol=0.01

    # Class fractions on the same disk
    M = [ismissing(d) || isnan(d) ? 0 : 1 for d in D]
    fr, fvm = remap_fractions(tgt, src, M, (0, 1); nsub=8)
    inside = fvm .== 1
    @test count(inside) > 1000
    @test fr[1][inside] ≈ fv[inside]
    @test all(fr[0][inside] .+ fr[1][inside] .≈ 1)

    # Same projection: sampled mean converges to the exact conservative mean
    fine = ProjGrid("F", (-400, 400), (-600, 200), 4, ps)
    c = ProjGrid("C", (-390, 390), (-590, 190), 20, ps)
    F = Float32.([sin(x / 50) * cos(y / 70) for x in fine.xc, y in fine.yc])
    Fe, _ = remap(c, fine, F)
    Fs, fvs = remap(c, fine, F; nsub=40)
    @test all(fvs .== 1) && maximum(abs.(Fs .- Fe)) < 0.02
end
