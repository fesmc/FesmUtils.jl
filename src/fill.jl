# Filling of gaps in fields: harmonic (Laplace) inpainting.
#
# The gap cells are the unknowns of the discrete Laplace equation (5-point stencil),
# with the valid cells next to them as fixed (Dirichlet) values. Neighbours that are
# neither gap nor valid (e.g. missing cells outside the gaps, or beyond the grid edge)
# are left out of the stencil (zero flux). The sparse system is solved directly.

"""
    laplace_fill!(F, gaps) -> nunfilled

Fill the cells of `F` where `gaps` is true by harmonic interpolation: the solution of
∇²u = 0 over the gaps (side neighbours, square cells), with the valid cells (not
NaN) next to them as boundary values; missing cells outside the gaps take no part.
Gaps not connected to any valid cell are set to NaN; returns their number of cells.
"""
function laplace_fill!(F::AbstractMatrix{<:AbstractFloat}, gaps::AbstractMatrix{Bool})
    size(gaps) == size(F) || throw(DimensionMismatch("field and gaps differ in size"))
    nx, ny = size(F)
    valid = map((v, g) -> !g && !isnan(v), F, gaps)
    inside(a, b) = 1 <= a <= nx && 1 <= b <= ny
    touches(i, j) = any(((di, dj),) -> inside(i + di, j + dj) && valid[i+di, j+dj], _NEIGHBOURS4)
    seeds = [gaps[i, j] && touches(i, j) for i in 1:nx, j in 1:ny]
    reached = flood_fill(gaps, seeds)
    unfilled = gaps .& .!reached
    F[unfilled] .= NaN
    cells = findall(reached)
    isempty(cells) && return count(unfilled)

    index = zeros(Int, nx, ny)
    for (k, I) in enumerate(cells)
        index[I] = k
    end
    rows, cols, vals = Int[], Int[], Float64[]
    b = zeros(length(cells))
    for (k, I) in enumerate(cells)
        i, j = Tuple(I)
        deg = 0
        for (di, dj) in _NEIGHBOURS4
            a, c = i + di, j + dj
            inside(a, c) || continue
            if reached[a, c]
                deg += 1
                push!(rows, k); push!(cols, index[a, c]); push!(vals, -1.0)
            elseif valid[a, c]
                deg += 1
                b[k] += F[a, c]
            end
        end
        push!(rows, k); push!(cols, k); push!(vals, deg)
    end
    u = sparse(rows, cols, vals, length(cells), length(cells)) \ b
    for (k, I) in enumerate(cells)
        F[I] = u[k]
    end
    return count(unfilled)
end

"""
    laplace_fill(F, gaps) -> (Ff, nunfilled)

`F` with its `gaps` filled by `laplace_fill!` (missing values as NaN), and the number
of gap cells not connected to any valid cell (NaN).
"""
function laplace_fill(F::AbstractMatrix, gaps::AbstractMatrix{Bool})
    T = _outtype(F)
    Ff = map(v -> _isvalid(v) ? T(v) : T(NaN), F)
    return Ff, laplace_fill!(Ff, gaps)
end
