# Rasterization of polygons onto a projected grid, by scanline filling at the cell
# centres: for each grid row, the crossings of the polygon edges with the row are
# sorted and the cells between pairs of crossings are filled (even-odd rule).

"""
    rasterize(g::ProjGrid, rings) -> BitMatrix

Cells of `g` whose centre lies inside the polygon given by `rings`, a vector of
rings of `(x, y)` vertices in the coordinates of `g` (km). Rings may be open or
closed. The even-odd rule is applied to all rings together, so the holes and parts
of a multipolygon are rings of the same list. Exact and threaded over rows.
"""
function rasterize(g::ProjGrid, rings::AbstractVector)
    nx, ny = size(g)
    x0, y0 = g.xc[1], g.yc[1]
    dx, dy = spacing(g)

    # Edges, and the rows whose centres they span (ymin <= y < ymax)
    edges = NTuple{4,Float64}[]
    for ring in rings
        n = length(ring)
        n < 3 && continue
        for k in 1:n
            (xa, ya) = ring[k]
            (xb, yb) = ring[k == n ? 1 : k + 1]
            ya == yb && continue
            push!(edges, (Float64(xa), Float64(ya), Float64(xb), Float64(yb)))
        end
    end
    rows = _edge_rows(edges, y0, dy, ny)

    # Edges crossing each row, in compressed form
    count = zeros(Int, ny + 1)
    for r in rows, j in r
        count[j+1] += 1
    end
    start = cumsum(count) .+ 1
    idx = Vector{Int}(undef, start[end] - 1)
    pos = start[1:ny]
    for (e, r) in enumerate(rows), j in r
        idx[pos[j]] = e
        pos[j] += 1
    end

    # Bool matrix while threaded: columns of a BitMatrix may share storage words
    mask = zeros(Bool, nx, ny)
    Threads.@threads for j in 1:ny
        y = y0 + (j - 1) * dy
        xs = Float64[]
        @inbounds for n in start[j]:start[j+1]-1
            xa, ya, xb, yb = edges[idx[n]]
            min(ya, yb) <= y < max(ya, yb) || continue
            push!(xs, xa + (y - ya) * (xb - xa) / (yb - ya))
        end
        sort!(xs)
        @inbounds for k in 1:2:length(xs)-1
            i0 = max(1, ceil(Int, (xs[k] - x0) / dx) + 1)
            i1 = min(nx, ceil(Int, (xs[k+1] - x0) / dx))
            for i in i0:i1
                mask[i, j] = true
            end
        end
    end
    return BitMatrix(mask)
end

# Rows j (centres y0 + (j-1) dy) that each edge may cross; checked exactly later.
function _edge_rows(edges, y0, dy, ny)
    return map(edges) do (xa, ya, xb, yb)
        lo, hi = minmax(ya, yb)
        j0 = max(1, ceil(Int, (lo - y0) / dy))
        j1 = min(ny, ceil(Int, (hi - y0) / dy) + 1)
        j0:j1
    end
end
