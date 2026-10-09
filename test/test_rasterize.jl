@testset "rasterize" begin
    proj = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +units=km"
    g = ProjGrid("t", (-10.0, 10.0), (-6.0, 8.0), 1.0, proj)

    # Brute force even-odd point-in-polygon at the cell centres
    function inside(x, y, rings)
        c = false
        for ring in rings, k in eachindex(ring)
            (xa, ya), (xb, yb) = ring[k], ring[k == length(ring) ? 1 : k + 1]
            if min(ya, yb) <= y < max(ya, yb) && x < xa + (y - ya) * (xb - xa) / (yb - ya)
                c = !c
            end
        end
        return c
    end
    brute(rings) = [inside(x, y, rings) for x in g.xc, y in g.yc]

    tri = [(-7.3, -4.2), (8.1, -1.7), (0.4, 6.6)]
    sq = [(-3.5, -3.5), (3.5, -3.5), (3.5, 3.5), (-3.5, 3.5), (-3.5, -3.5)]   # closed
    hole = [(-1.5, -1.5), (1.5, -1.5), (1.5, 1.5), (-1.5, 1.5)]
    far = [(20.0, 20.0), (30.0, 20.0), (25.0, 30.0)]
    for rings in ([tri], [sq], [sq, hole], [tri, far], [sq, [(5.2, 2.1), (9.7, 2.1), (9.7, 7.3)]])
        @test rasterize(g, rings) == brute(rings)
    end
    # Cells centred on a horizontal edge: lower edge inside, upper outside
    m = rasterize(g, [[(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)]])
    @test count(m) == 16
    @test !any(rasterize(g, [far]))
end
