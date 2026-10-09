@testset "remap aligned" begin
    proj = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    fine = ProjGrid("F", (-400, 400), (-600, 200), 4, proj)
    F = rand(Float32, size(fine)...)

    # Identity
    Ft, fv = remap(fine, fine, F)
    @test Ft ≈ F
    @test all(fv .≈ 1)

    # 2x node-nested coarsening = 1/4-1/2-1/4 weighting in the interior
    c = coarsen(fine, 2; name="C")
    Ft, fv = remap(c, fine, F)
    w = [0.25, 0.5, 0.25]
    i, j = 10, 20
    ref = sum(w[a] * w[b] * F[2i-3+a, 2j-3+b] for a in 1:3, b in 1:3)
    @test Ft[i, j] ≈ ref
    @test fv[1, 1] ≈ 0.75^2 && fv[2, 2] ≈ 1   # corner cell extends dx/4 beyond the fine grid on each side

    # Conservation onto a target that covers the source, with non-integer ratio
    big = ProjGrid("B", (-410, 410), (-610, 210), 10, proj)
    Ft, fv = remap(big, fine, F)
    @test sum(Float64.(Ft[fv .> 0]) .* fv[fv .> 0]) * 100 ≈ sum(Float64.(F)) * 16 rtol=1e-6

    # Missing values are excluded
    Fm = Array{Union{Missing,Float32}}(F)
    Fm[1:2, 1:2] .= missing
    Ft, fv = remap(c, fine, Fm)
    @test fv[1, 1] == 0 && isnan(Ft[1, 1])

    # Threads give the same answer as serial weights applied directly
    m = AlignedMap(c, fine)
    @test remap(m, F)[1] == remap(c, fine, F)[1]

    # Different projection is rejected
    south = ProjGrid("S", (-400, 400), (-600, 200), 4,
                     polar_stereographic_proj(lat_0=-90, lat_ts=-71, lon_0=0, a=6378137, rf=298.257223563))
    @test_throws ArgumentError remap(south, fine, F)

    # Class fractions sum to one
    M = rand(1:3, size(fine)...)
    fr, fv = remap_fractions(c, fine, M, 1:3)
    @test all(fr[1] .+ fr[2] .+ fr[3] .≈ 1)
end

@testset "remap_dominant" begin
    proj = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km"
    src = ProjGrid("s", (0.0, 11.0), (0.0, 7.0), 1.0, proj)
    M = [i <= 5 ? 1 : (j <= 3 ? 2 : 3) for i in 1:12, j in 1:8]
    # Factor 1: unchanged
    @test remap_dominant(src, src, M) == M
    # Odd factor: coarse cells are 3x3 blocks centred on fine cells
    tgt = coarsen(src, 3; name="t")
    Mt = remap_dominant(tgt, src, M)
    for i in axes(Mt, 1), j in axes(Mt, 2)
        ic, jc = 3(i - 1) + 1, 3(j - 1) + 1
        blk = [M[k, l] for k in max(1, ic - 1):min(12, ic + 1), l in max(1, jc - 1):min(8, jc + 1)]
        cnt = Dict(c => count(==(c), blk) for c in unique(blk))
        best = maximum(values(cnt))
        @test Mt[i, j] == minimum(c for (c, n) in cnt if n == best)
    end
end
