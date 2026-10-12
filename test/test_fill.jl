@testset "laplace_fill" begin
    # A linear field is harmonic: a hole in it is filled exactly
    nx, ny = 30, 20
    L = Float32[2i - 3j + 1 for i in 1:nx, j in 1:ny]
    gaps = falses(nx, ny)
    gaps[10:15, 5:12] .= true
    F = copy(L)
    F[gaps] .= NaN32
    @test laplace_fill!(F, gaps) == 0
    @test F ≈ L
    Ff, n = laplace_fill(L .+ 100 .* gaps, gaps)       # values in the gaps are replaced
    @test n == 0 && Ff ≈ L

    # Missing cells outside the gaps take no part (zero flux): a gap against a missing
    # region in a constant field stays constant; the maximum principle holds
    F = fill(5.0, nx, ny)
    F[1:12, :] .= NaN
    g = falses(nx, ny)
    g[13:16, 3:6] .= true
    F[g] .= NaN
    @test laplace_fill!(F, g) == 0
    @test all(F[g] .≈ 5)
    R = rand(nx, ny)
    R[g] .= NaN
    laplace_fill!(R, g)
    @test all(0 .<= R[g] .<= 1)

    # Gaps not connected to valid cells stay missing and are counted; `missing` input
    M = Union{Missing,Float32}[1.0f0 for i in 1:nx, j in 1:ny]
    M[1:5, :] .= missing
    g = falses(nx, ny)
    g[2:3, 2:3] .= true                         # surrounded by missing cells
    g[20, 10] = true                            # surrounded by valid cells
    Mf, n = laplace_fill(M, g)
    @test n == 4
    @test eltype(Mf) == Float32
    @test all(isnan, Mf[2:3, 2:3]) && Mf[20, 10] ≈ 1
    @test all(isnan, Mf[1:5, 1]) && Mf[6, 1] == 1

    # A gap at the grid edge (zero flux there), and the mean of the neighbours for a
    # single cell
    F = Float64[i for i in 1:5, j in 1:5]
    g = falses(5, 5)
    g[1, 1:5] .= true
    F[g] .= NaN
    laplace_fill!(F, g)
    @test all(F[1, :] .≈ 2)
    S = [0.0 1 0; 2 NaN 4; 0 8 0]
    s = falses(3, 3)
    s[2, 2] = true
    laplace_fill!(S, s)
    @test S[2, 2] ≈ (1 + 2 + 4 + 8) / 4
    @test_throws DimensionMismatch laplace_fill!(zeros(3, 3), falses(2, 3))
end
