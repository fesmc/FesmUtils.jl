# Exact Euclidean distance transform on a regular grid (Felzenszwalb & Huttenlocher,
# 2012, Theory of Computing 8:415-428): two separable passes of the 1D transform.

"""
    distance_to(mask, dx, dy=dx) -> D

Euclidean distance from each cell centre to the nearest cell centre where `mask` is
true, for cell spacings `dx` and `dy` (0 inside `mask`, Inf if `mask` is all
false). Exact and threaded.
"""
function distance_to(mask::AbstractMatrix{Bool}, dx::Real, dy::Real=dx)
    nx, ny = size(mask)
    G = Matrix{Float64}(undef, nx, ny)
    Threads.@threads for j in 1:ny
        f = [mask[i, j] ? 0.0 : Inf for i in 1:nx]
        G[:, j] = _edt1d(f, Float64(dx))
    end
    D = Matrix{Float64}(undef, nx, ny)
    Threads.@threads for i in 1:nx
        D[i, :] = sqrt.(_edt1d(G[i, :], Float64(dy)))
    end
    return D
end

# Squared distance transform of sampled function f (Inf = no point) with spacing h:
# d(q) = min_p (h (q - p))^2 + f(p), from the lower envelope of parabolas.
function _edt1d(f::Vector{Float64}, h::Float64)
    n = length(f)
    h2 = h^2
    v = zeros(Int, n)          # parabola locations in the envelope
    z = zeros(Float64, n + 1)  # boundaries between parabolas
    k = 0
    for q in 1:n
        isfinite(f[q]) || continue
        while k > 0
            p = v[k]
            s = ((f[q] + h2 * q^2) - (f[p] + h2 * p^2)) / (2h2 * (q - p))
            if s <= z[k]
                k -= 1
            else
                k += 1
                v[k] = q
                z[k] = s
                z[k+1] = Inf
                break
            end
        end
        if k == 0
            k = 1
            v[1] = q
            z[1] = -Inf
            z[2] = Inf
        end
    end

    d = fill(Inf, n)
    k == 0 && return d
    k = 1
    for q in 1:n
        while z[k+1] < q
            k += 1
        end
        d[q] = h2 * (q - v[k])^2 + f[v[k]]
    end
    return d
end
