# Bilinear interpolation between grids: each target cell takes the value at its
# centre, interpolated between the four surrounding source cell centres.

_sample_transform(::LonLatGrid, ::LonLatGrid) = identity

_centres(g::ProjGrid) = (g.xc, g.yc)
_centres(g::LonLatGrid) = (g.lon, g.lat)

# Source cells i0, i1 = i0 + 1 (periodically) around the position f (in cells from
# the first centre, 0-based) along an axis of n cells, and the weight t of i1; nothing
# outside the axis. Between the outermost centres and the edge, the outermost cell
# is used.
function _bracket(f::Float64, n::Int, periodic::Bool)
    if periodic
        i0 = floor(Int, f)
        return (mod(i0, n) + 1, mod(i0 + 1, n) + 1, f - i0)
    end
    -0.5 <= f <= n - 0.5 || return nothing
    f = clamp(f, 0.0, n - 1.0)
    i0 = min(floor(Int, f), n - 2)
    return (i0 + 1, i0 + 2, f - i0)
end

_position(src::ProjGrid, x, y) =
    (isfinite(x) && isfinite(y)) ? (_bracket((x - src.xc[1]) / spacing(src)[1], length(src.xc), false),
                                    _bracket((y - src.yc[1]) / spacing(src)[2], length(src.yc), false)) :
                                   (nothing, nothing)

function _position(src::LonLatGrid, lon, lat)
    (isfinite(lon) && isfinite(lat)) || return (nothing, nothing)
    n = length(src.lon)
    fx = src.isglobal ? mod(lon - src.lon[1], 360.0) / src.dlon :
                        mod(lon - src.lon[1] + src.dlon / 2, 360.0) / src.dlon - 0.5
    return (_bracket(fx, n, src.isglobal), _bracket((lat - src.lat[1]) / src.dlat, length(src.lat), false))
end

"""
    remap_bilinear(tgt, src, F) -> (Ft, f_valid)

Bilinear interpolation of the field `F` (size of `src`) at the cell centres of `tgt`.
The grids are projected (`ProjGrid`, any projection) or lon-lat (`LonLatGrid`, periodic
in longitude when global). Missing source values (`missing` or NaN) are left out and
the weights of the others renormalised; `f_valid` is the weight of the valid source
cells, 0 outside the source (NaN in `Ft`). Threaded over target rows.
"""
function remap_bilinear(tgt::Union{ProjGrid,LonLatGrid}, src::Union{ProjGrid,LonLatGrid}, F::AbstractMatrix)
    size(F) == size(src) ||
        throw(DimensionMismatch("field size $(size(F)) does not match source grid $(size(src))"))
    T = _outtype(F)
    xs, ys = _centres(tgt)
    Ft = fill(T(NaN), size(tgt))
    fv = zeros(Float32, size(tgt))
    tasks = map(_chunks(length(ys))) do js
        Threads.@spawn begin
            trans = _sample_transform(tgt, src)
            @inbounds for j in js, i in eachindex(xs)
                u, v = trans((xs[i], ys[j]))
                bx, by = _position(src, u, v)
                (bx === nothing || by === nothing) && continue
                (i0, i1, tx), (j0, j1, ty) = bx, by
                a = 0.0
                b = 0.0
                for (l, k, w) in ((i0, j0, (1 - tx) * (1 - ty)), (i1, j0, tx * (1 - ty)),
                                  (i0, j1, (1 - tx) * ty), (i1, j1, tx * ty))
                    val = F[l, k]
                    if _isvalid(val) && w > 0
                        a += w * val
                        b += w
                    end
                end
                if b > 0
                    Ft[i, j] = T(a / b)
                    fv[i, j] = Float32(b)
                end
            end
        end
    end
    foreach(wait, tasks)
    return Ft, fv
end
