@testset "smooth" begin
    proj = polar_stereographic_proj(lat_0=90, lat_ts=70, lon_0=-45, a=6378137, rf=298.257223563)
    g = ProjGrid("G", (-200, 200), (-300, 100), 4, proj)
    nx, ny = size(g)

    # A constant stays constant, also next to missing cells; no smoothing is the identity
    F = fill(3.0f0, nx, ny)
    F[40:50, 40:50] .= NaN32
    Fs = smooth(g, F, 20)
    @test eltype(Fs) == Float32
    @test all(isnan, Fs[40:50, 40:50])
    @test all(Fs[.!isnan.(F)] .≈ 3)
    R = rand(Float32, nx, ny)
    @test smooth(g, R, 0) ≈ R

    # A point source spreads as a Gaussian with the given standard deviation (km)
    D = zeros(nx, ny)
    D[51, 51] = 1
    Ds = smooth(g, D, 12)
    @test sum(Ds) ≈ 1
    xs = (1:nx) .- 51
    @test sqrt(sum(Ds[:, 51] .* (4xs) .^ 2) / sum(Ds[:, 51])) ≈ 12 rtol=0.02
    @test Ds[51 + 3, 51] ≈ Ds[51, 51 + 3]

    # Lon-lat: periodic in longitude on a global grid, constant preserved
    ll = LonLatGrid(-179.5:1:179.5, -89.5:1:89.5)
    L = zeros(size(ll))
    L[1, 90] = 1
    Ls = smooth(ll, L, 300)
    @test Ls[end, 90] ≈ Ls[2, 90]
    @test Ls[end, 90] > 0
    @test all(smooth(ll, fill(2.0, size(ll)), 500) .≈ 2)

    # Lon-lat kernel is wider in degrees of longitude at high latitude
    P = zeros(size(ll))
    P[180, 10] = 1                                   # 80.5°S
    P[180, 90] = 1                                   # 0.5°S
    Ps = smooth(ll, P, 200)
    @test count(>(1e-3 * maximum(Ps[:, 10])), Ps[:, 10]) > 3 * count(>(1e-3 * maximum(Ps[:, 90])), Ps[:, 90])
end
