# Remapping from a regular lon-lat grid onto a projected grid by supersampling.
#
# Each target cell is sampled at nsub x nsub points spaced uniformly in the
# projected plane. Each sample takes the value of the source cell that contains it,
# so the mean over samples converges to the conservative (area-weighted) cell mean
# as nsub grows. Choose nsub so the sample spacing is at most about half the
# source spacing. Threaded over target rows, with one Proj transformation per task.

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

# Call f(i, j, is, js) for every sample of every target cell (i, j), where
# (is, js) is the source cell under the sample. Calls for a given (i, j) all happen
# on the same task, so f may update target arrays at (i, j) without locking.
function _supersample(f, tgt::ProjGrid, src::LonLatGrid, nsub::Integer)
    nsub >= 1 || throw(ArgumentError("nsub must be >= 1"))
    nx, ny = size(tgt)
    dx, dy = spacing(tgt)
    off = ((1:nsub) .- (nsub + 1) / 2) ./ nsub
    tasks = map(_chunks(ny)) do js
        Threads.@spawn begin
            trans = _inverse_transform(tgt)
            for j in js, i in 1:nx, oy in off, ox in off
                lon, lat = trans((tgt.xc[i] + ox * dx, tgt.yc[j] + oy * dy))
                is, jsrc = _cell_index(src, lon, lat)
                f(i, j, is, jsrc)
            end
        end
    end
    foreach(wait, tasks)
    return nothing
end

"""
    remap(tgt::ProjGrid, src::LonLatGrid, F; nsub) -> (Ft, f_valid)

Area-weighted mean of the lon-lat field `F` (size of `src`) over each target cell,
estimated from `nsub` x `nsub` samples per cell. `Ft` is the mean over valid samples
(NaN where there are none), and `f_valid` the fraction of valid samples. Threaded.
"""
function remap(tgt::ProjGrid, src::LonLatGrid, F::AbstractMatrix; nsub::Integer)
    size(F) == size(src) ||
        throw(DimensionMismatch("field size $(size(F)) does not match source grid $(size(src))"))
    num = zeros(Float64, size(tgt))
    cnt = zeros(Int32, size(tgt))
    _supersample(tgt, src, nsub) do i, j, is, js
        is == 0 && return
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

Area fraction of each class of the categorical lon-lat field `M` within each target
cell, relative to its valid part, from `nsub` x `nsub` samples per cell.
"""
function remap_fractions(tgt::ProjGrid, src::LonLatGrid, M::AbstractMatrix, classes; nsub::Integer)
    size(M) == size(src) ||
        throw(DimensionMismatch("field size $(size(M)) does not match source grid $(size(src))"))
    cls = collect(classes)
    cnt = zeros(Int32, length(cls), size(tgt)...)
    nvalid = zeros(Int32, size(tgt))
    _supersample(tgt, src, nsub) do i, j, is, js
        is == 0 && return
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
