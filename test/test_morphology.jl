@testset "morphology" begin
    m = falses(21, 21)
    m[6:16, 6:16] .= true
    @test erode(m, 2.0, 1.0) == (falses(21, 21) .| [(8 <= i <= 14 && 8 <= j <= 14) for i in 1:21, j in 1:21])
    d = dilate(m, 1.0, 1.0)
    @test count(d) == 11 * 11 + 4 * 11
    @test dilate(erode(m, 2.0, 1.0), 2.0, 1.0) ⊆ m

    # Flood fill: two blobs, one seeded
    m = falses(10, 10)
    m[1:3, 1:3] .= true
    m[6:9, 6:9] .= true
    s = falses(10, 10)
    s[2, 2] = true
    f = flood_fill(m, s)
    @test count(f) == 9 && f[1:3, 1:3] == trues(3, 3)
    @test !any(flood_fill(m, .!m))

    # Label extension: through allowed cells only, nearest along paths
    L = zeros(Int, 9, 5)
    L[1, 3] = 1
    L[9, 3] = 2
    allowed = trues(9, 5)
    E = extend_labels(L, allowed, 1.0)
    @test all(E[1:4, :] .== 1) && all(E[6:9, :] .== 2)
    # A wall at column 6 except one gap far away: cells right of the wall near label 1
    # are reached only around it
    allowed[6, 1:4] .= false
    L[9, 3] = 0
    L[9, 1] = 2
    E = extend_labels(L, allowed, 1.0)
    @test E[6, 1] == 0 && E[7, 3] == 2
    @test extend_labels(zeros(Int, 3, 3), trues(3, 3), 1.0) == zeros(Int, 3, 3)
end
