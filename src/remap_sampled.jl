# Remapping onto a projected grid by supersampling, from a regular lon-lat grid or
# from a projected grid on another projection.
#
# Each target cell is sampled at nsub x nsub points spaced uniformly in the
# projected plane. Each sample takes the value of the source cell that contains it,
# so the mean over samples converges to the conservative (area-weighted) cell mean
# as nsub grows. Choose nsub so the sample spacing is at most about half the
# source spacing. Threaded over target rows, with one Proj transformation per task.

# Transformation from target to source coordinates, created per task (Proj
# transformations are not thread-safe).
_sample_transform(tgt::ProjGrid, ::LonLatGrid) = _inverse_transform(tgt)
_sample_transform(tgt::ProjGrid, src::ProjGrid) =
    Proj.Transformation(tgt.proj, src.proj; always_xy=true, ctx=Proj.proj_context_create())

# Source cell containing (lon, lat), or (0, 0) if outside the source grid.
@inline function _cell_index(src::LonLatGrid, lon::Float64, lat::Float64)
    j = floor(Int, (lat - src.lat[1]) / src.dlat + 0.5) + 1
    1 <= j <= length(src.lat) || return (0, 0)
    i = floor(Int, mod(lon - src.lon[1] + src.dlon / 2, 360.0) / src.dlon) + 1
    if i > length(src.lon)
        src.isglobal || return (0, 0)
        i = 1
    end
    return (i, j)
end

# Source cell containing (x, y), or (0, 0) if outside the source grid or not
# representable in its projection.
@inline function _cell_index(src::ProjGrid, x::Float64, y::Float64)
    (isfinite(x) && isfinite(y)) || return (0, 0)
    dx, dy = spacing(src)
    fi = (x - src.xc[1]) / dx + 0.5
    fj = (y - src.yc[1]) / dy + 0.5
    (0 <= fi < length(src.xc) && 0 <= fj < length(src.yc)) || return (0, 0)
    return (floor(Int, fi) + 1, floor(Int, fj) + 1)
end

# Outline of a projected source grid in target coordinates: a closed polygon
# (px, py), and the distance m by which it is widened so that it contains every
# target point inside the source. The outline is curved on the target, so m is
# twice the largest offset of the outline from the polygon edges (at their
# midpoints), plus one sample spacing to cover rounding and the inexact inverse
# transformation. `nothing` for lon-lat sources, or if the outline is not
# representable on the target projection.
_footprint(tgt::ProjGrid, src::LonLatGrid, nsub::Integer) = nothing

function _footprint(tgt::ProjGrid, src::ProjGrid, nsub::Integer; npts::Integer=32)
    dx, dy = spacing(src)
    x0, x1 = src.xc[1] - dx / 2, src.xc[end] + dx / 2
    y0, y1 = src.yc[1] - dy / 2, src.yc[end] + dy / 2
    # Vertices of the polygon, alternating with the midpoints of its edges
    t = (0:2npts-1) ./ (2npts)
    pts = vcat([(x0 + s * (x1 - x0), y0) for s in t], [(x1, y0 + s * (y1 - y0)) for s in t],
               [(x1 - s * (x1 - x0), y1) for s in t], [(x0, y1 - s * (y1 - y0)) for s in t])
    trans = Proj.Transformation(src.proj, tgt.proj; always_xy=true, ctx=Proj.proj_context_create())
    xy = map(trans, pts)
    all(p -> isfinite(p[1]) && isfinite(p[2]), xy) || return nothing
    push!(xy, xy[1])
    sag = 0.0
    for k in 1:2:length(xy)-2
        (ax, ay), (mx, my), (bx, by) = xy[k], xy[k+1], xy[k+2]
        sag = max(sag, abs((bx - ax) * (my - ay) - (by - ay) * (mx - ax)) / hypot(bx - ax, by - ay))
    end
    tdx, tdy = spacing(tgt)
    return first.(xy[1:2:end]), last.(xy[1:2:end]), 2sag + max(tdx, tdy) / nsub
end

# Intervals of x, sorted and disjoint, where the line y = yr is inside the polygon
# (px, py) or within m of an edge.
function _row_intervals!(ivs::Vector{NTuple{2,Float64}}, cross::Vector{Float64}, fp, yr::Float64)
    px, py, m = fp
    empty!(ivs)
    empty!(cross)
    for k in 1:length(px)-1
        ax, ay, bx, by = px[k], py[k], px[k+1], py[k+1]
        # Crossing of the line (half-open in y, so each crossing counts once)
        (ay <= yr) != (by <= yr) && push!(cross, ax + (yr - ay) * (bx - ax) / (by - ay))
        # Part of the edge within m of the line in y, widened by m in x
        if ay == by
            abs(yr - ay) <= m && push!(ivs, (min(ax, bx) - m, max(ax, bx) + m))
        else
            t1, t2 = minmax((yr - m - ay) / (by - ay), (yr + m - ay) / (by - ay))
            t1, t2 = max(t1, 0.0), min(t2, 1.0)
            if t1 <= t2
                xa, xb = minmax(ax + t1 * (bx - ax), ax + t2 * (bx - ax))
                push!(ivs, (xa - m, xb + m))
            end
        end
    end
    sort!(cross)
    for k in 1:2:length(cross)-1
        push!(ivs, (cross[k], cross[k+1]))
    end
    sort!(ivs)
    n = 0
    for iv in ivs
        if n > 0 && iv[1] <= ivs[n][2]
            ivs[n] = (ivs[n][1], max(ivs[n][2], iv[2]))
        else
            n += 1
            ivs[n] = iv
        end
    end
    resize!(ivs, n)
    return ivs
end

# Call f(i, j, is, js) for every sample of every target cell (i, j) that falls in
# a source cell (is, js). Calls for a given (i, j) all happen on the same task, in
# the same order (rows of samples, then along each row), so f may update target
# arrays at (i, j) without locking, and sums are reproducible. For a projected
# source, samples outside its outline on the target are skipped before they are
# transformed: the result is the same, but a small source on a large target needs
# far fewer transformations.
function _supersample(f, tgt::ProjGrid, src::Union{LonLatGrid,ProjGrid}, nsub::Integer)
    nsub >= 1 || throw(ArgumentError("nsub must be >= 1"))
    nx, ny = size(tgt)
    dx, dy = spacing(tgt)
    off = ((1:nsub) .- (nsub + 1) / 2) ./ nsub
    fp = _footprint(tgt, src, nsub)
    tasks = map(_chunks(ny)) do js
        Threads.@spawn begin
            trans = _sample_transform(tgt, src)
            ivs = [(-Inf, Inf)]
            cross = Float64[]
            for j in js, oy in off
                y = tgt.yc[j] + oy * dy
                fp === nothing || _row_intervals!(ivs, cross, fp, y)
                isempty(ivs) && continue
                p = 1
                for i in 1:nx, ox in off
                    x = tgt.xc[i] + ox * dx
                    while p <= length(ivs) && x > ivs[p][2]
                        p += 1
                    end
                    p > length(ivs) && break
                    x < ivs[p][1] && continue
                    u, v = trans((x, y))
                    is, jsrc = _cell_index(src, u, v)
                    is == 0 || f(i, j, is, jsrc)
                end
            end
        end
    end
    foreach(wait, tasks)
    return nothing
end

"""
    remap(tgt::ProjGrid, src::LonLatGrid, F; nsub) -> (Ft, f_valid)
    remap(tgt::ProjGrid, src::ProjGrid, F; nsub) -> (Ft, f_valid)

Area-weighted mean of the field `F` (size of `src`) over each target cell,
estimated from `nsub` x `nsub` samples per cell. The source is a lon-lat grid, or a
projected grid on any projection. `Ft` is the mean over valid samples (NaN where
there are none), and `f_valid` the fraction of valid samples. Threaded.
"""
remap(tgt::ProjGrid, src::LonLatGrid, F::AbstractMatrix; nsub::Integer) =
    _remap_sampled(tgt, src, F, nsub)

function _remap_sampled(tgt::ProjGrid, src, F::AbstractMatrix, nsub::Integer)
    size(F) == size(src) ||
        throw(DimensionMismatch("field size $(size(F)) does not match source grid $(size(src))"))
    num = zeros(Float64, size(tgt))
    cnt = zeros(Int32, size(tgt))
    _supersample(tgt, src, nsub) do i, j, is, js
        @inbounds v = F[is, js]
        _isvalid(v) || return
        @inbounds num[i, j] += v
        @inbounds cnt[i, j] += 1
    end
    T = _outtype(F)
    Ft = map((a, n) -> n > 0 ? T(a / n) : T(NaN), num, cnt)
    fv = Float32.(cnt ./ nsub^2)
    return Ft, fv
end

"""
    remap_fractions(tgt::ProjGrid, src::LonLatGrid, M, classes; nsub) -> (fracs, f_valid)
    remap_fractions(tgt::ProjGrid, src::ProjGrid, M, classes; nsub) -> (fracs, f_valid)

Area fraction of each class of the categorical field `M` (on a lon-lat grid, or a
projected grid on any projection) within each target cell, relative to its valid
part, from `nsub` x `nsub` samples per cell.
"""
remap_fractions(tgt::ProjGrid, src::LonLatGrid, M::AbstractMatrix, classes; nsub::Integer) =
    _remap_fractions_sampled(tgt, src, M, classes, nsub)

function _remap_fractions_sampled(tgt::ProjGrid, src, M::AbstractMatrix, classes, nsub::Integer)
    size(M) == size(src) ||
        throw(DimensionMismatch("field size $(size(M)) does not match source grid $(size(src))"))
    cls = collect(classes)
    cnt = zeros(Int32, length(cls), size(tgt)...)
    nvalid = zeros(Int32, size(tgt))
    _supersample(tgt, src, nsub) do i, j, is, js
        @inbounds v = M[is, js]
        ismissing(v) && return
        @inbounds nvalid[i, j] += 1
        c = findfirst(==(v), cls)
        c === nothing || (@inbounds cnt[c, i, j] += 1)
    end
    fracs = Dict{eltype(cls),Matrix{Float32}}()
    for (c, cl) in enumerate(cls)
        fracs[cl] = map((n, nv) -> nv > 0 ? Float32(n / nv) : NaN32, view(cnt, c, :, :), nvalid)
    end
    return fracs, Float32.(nvalid ./ nsub^2)
end
