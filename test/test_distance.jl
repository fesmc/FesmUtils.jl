@testset "distance" begin
    # Brute force reference
    function brute(mask, dx, dy)
        pts = [(i, j) for i in axes(mask, 1), j in axes(mask, 2) if mask[i, j]]
        isempty(pts) && return fill(Inf, size(mask))
        return [minimum(sqrt(((i - p) * dx)^2 + ((j - q) * dy)^2) for (p, q) in pts)
                for i in axes(mask, 1), j in axes(mask, 2)]
    end
    for (nx, ny, frac) in ((37, 23, 0.02), (50, 50, 0.2), (8, 30, 0.5))
        mask = rand(nx, ny) .< frac
        mask[1, 1] = true
        @test distance_to(mask, 2.0, 3.0) ≈ brute(mask, 2.0, 3.0)
    end
    @test all(isinf, distance_to(falses(5, 4), 1.0))
    @test distance_to(trues(5, 4), 1.0) == zeros(5, 4)
end
