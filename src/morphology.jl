# Operations on masks and label fields on regular grids: erosion and dilation by a
# distance, flood filling from seed cells, and extension of labels along paths.

"""
    erode(mask, r, dx, dy=dx) -> BitMatrix

Cells of `mask` farther than `r` from every cell outside it (same units as `dx`).
Cells beyond the grid edge do not count as outside.
"""
erode(mask::AbstractMatrix{Bool}, r::Real, dx::Real, dy::Real=dx) =
    distance_to(.!mask, dx, dy) .> r

"""
    dilate(mask, r, dx, dy=dx) -> BitMatrix

Cells within distance `r` of a cell of `mask` (same units as `dx`).
"""
dilate(mask::AbstractMatrix{Bool}, r::Real, dx::Real, dy::Real=dx) =
    distance_to(mask, dx, dy) .<= r

const _NEIGHBOURS4 = ((1, 0), (-1, 0), (0, 1), (0, -1))

"""
    flood_fill(mask, seeds) -> BitMatrix

Cells of `mask` connected to a cell of `mask .& seeds` through cells of `mask`,
with side neighbours (4-connectivity).
"""
function flood_fill(mask::AbstractMatrix{Bool}, seeds::AbstractMatrix{Bool})
    size(seeds) == size(mask) || throw(DimensionMismatch("mask and seeds differ in size"))
    nx, ny = size(mask)
    filled = mask .& seeds
    stack = [(i, j) for i in 1:nx, j in 1:ny if filled[i, j]]
    @inbounds while !isempty(stack)
        i, j = pop!(stack)
        for (di, dj) in _NEIGHBOURS4
            a, b = i + di, j + dj
            (1 <= a <= nx && 1 <= b <= ny) || continue
            if mask[a, b] && !filled[a, b]
                filled[a, b] = true
                push!(stack, (a, b))
            end
        end
    end
    return filled
end

const _NEIGHBOURS8 = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))

"""
    extend_labels(L, allowed, dx, dy=dx) -> L2

Extend the labels `L` (0 = no label) into the unlabelled cells of `allowed`: each
such cell gets the label of the nearest labelled cell, with distances measured
along paths through `allowed` cells (steps to the 8 neighbours). Cells that no path
reaches keep 0.
"""
function extend_labels(L::AbstractMatrix{T}, allowed::AbstractMatrix{Bool},
                       dx::Real, dy::Real=dx) where {T<:Integer}
    size(allowed) == size(L) || throw(DimensionMismatch("labels and allowed differ in size"))
    nx, ny = size(L)
    L2 = Matrix{T}(L)
    D = fill(Inf, nx, ny)
    steps = map(((di, dj),) -> hypot(di * dx, dj * dy), _NEIGHBOURS8)

    # Dijkstra from the labelled cells next to cells to fill, with a binary heap of
    # (distance, i, j)
    fillable(a, b) = 1 <= a <= nx && 1 <= b <= ny && allowed[a, b] && L[a, b] == 0
    heap = Tuple{Float64,Int,Int}[]
    for j in 1:ny, i in 1:nx
        L[i, j] == 0 && continue
        D[i, j] = 0.0
        any(((di, dj),) -> fillable(i + di, j + dj), _NEIGHBOURS8) && _heap_push!(heap, (0.0, i, j))
    end
    @inbounds while !isempty(heap)
        d, i, j = _heap_pop!(heap)
        d > D[i, j] && continue
        for (n, (di, dj)) in enumerate(_NEIGHBOURS8)
            a, b = i + di, j + dj
            fillable(a, b) || continue
            dn = d + steps[n]
            if dn < D[a, b]
                D[a, b] = dn
                L2[a, b] = L2[i, j]
                _heap_push!(heap, (dn, a, b))
            end
        end
    end
    return L2
end

function _heap_push!(h::Vector, v)
    push!(h, v)
    k = length(h)
    @inbounds while k > 1
        p = k >> 1
        h[p] <= h[k] && break
        h[p], h[k] = h[k], h[p]
        k = p
    end
    return h
end

function _heap_pop!(h::Vector)
    top = h[1]
    last = pop!(h)
    n = length(h)
    if n > 0
        h[1] = last
        k = 1
        @inbounds while true
            c = 2k
            c > n && break
            (c < n && h[c+1] < h[c]) && (c += 1)
            h[k] <= h[c] && break
            h[k], h[c] = h[c], h[k]
            k = c
        end
    end
    return top
end
