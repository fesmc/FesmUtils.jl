# Exact conservative remapping between regular grids on the same projection.
#
# Cells of both grids are axis-aligned rectangles, so the overlap area of a source
# and a target cell is the product of their 1D overlaps along x and y. The weights
# are then two small 1D operators, applied in two threaded passes:
#
#     Ft = Wx * F * Wy'
#
# Missing source values (`missing` or NaN) are excluded, and the result is the mean
# over the valid part of each target cell. Areas are projected areas.

"""
    Overlap1D(tgt, src)

Overlap weights of 1D cells with centres `tgt` onto cells with centres `src`:
for target cell `i`, source cells `ranges[i]` cover the fractions `weights[i]`
of its length.
"""
struct Overlap1D
    ranges::Vector{UnitRange{Int}}
    weights::Vector{Vector{Float64}}
end

function Overlap1D(tgt::AbstractVector{<:Real}, src::AbstractVector{<:Real})
    dt = (tgt[end] - tgt[1]) / (length(tgt) - 1)
    ds = (src[end] - src[1]) / (length(src) - 1)
    ns = length(src)
    tol = 1e-9 * min(dt, ds)

    ranges = Vector{UnitRange{Int}}(undef, length(tgt))
    weights = Vector{Vector{Float64}}(undef, length(tgt))
    for (i, xt) in enumerate(tgt)
        lo, hi = xt - dt / 2, xt + dt / 2
        k0 = max(1, floor(Int, (lo - (src[1] - ds / 2)) / ds) + 1)
        k1 = min(ns, ceil(Int, (hi - (src[1] - ds / 2)) / ds))
        ks = Int[]
        ws = Float64[]
        for k in k0:k1
            ov = min(hi, src[k] + ds / 2) - max(lo, src[k] - ds / 2)
            if ov > tol
                push!(ks, k)
                push!(ws, ov / dt)
            end
        end
        ranges[i] = isempty(ks) ? (1:0) : (ks[1]:ks[end])
        weights[i] = ws
    end
    return Overlap1D(ranges, weights)
end

"""
    AlignedMap(tgt, src)

Conservative remapping weights from `src` to `tgt`, two `ProjGrid`s on the same
projection. Build once and reuse for all fields of a source.
"""
struct AlignedMap
    tgt::ProjGrid
    src::ProjGrid
    wx::Overlap1D
    wy::Overlap1D
end

function AlignedMap(tgt::ProjGrid, src::ProjGrid)
    same_projection(tgt, src) ||
        throw(ArgumentError("grids $(tgt.name) and $(src.name) are not on the same projection"))
    return AlignedMap(tgt, src, Overlap1D(tgt.xc, src.xc), Overlap1D(tgt.yc, src.yc))
end

"""
    same_projection(a, b)

True if `a` and `b` use the same map projection, tested by transforming points
across the domain of `b` and checking that they are unchanged.
"""
function same_projection(a::ProjGrid, b::ProjGrid)
    a.proj == b.proj && return true
    trans = Proj.Transformation(b.proj, a.proj; always_xy=true)
    for x in (b.xc[1], b.xc[end], 0.5 * (b.xc[1] + b.xc[end])),
        y in (b.yc[1], b.yc[end], 0.5 * (b.yc[1] + b.yc[end]))
        xa, ya = trans((x, y))
        (isapprox(xa, x; atol=1e-6) && isapprox(ya, y; atol=1e-6)) || return false
    end
    return true
end

_isvalid(v) = !ismissing(v) && !isnan(v)

# Output element type: keep the input float precision, Float32 otherwise.
function _outtype(F::AbstractArray)
    T = nonmissingtype(eltype(F))
    return T <: AbstractFloat ? T : Float32
end

"""
    remap(tgt::ProjGrid, src::ProjGrid, F) -> (Ft, f_valid)
    remap(m::AlignedMap, F) -> (Ft, f_valid)

Exact conservative remapping of `F` (size of `src`) onto `tgt`, on the same
projection. `Ft` is the mean over the valid part of each target cell (NaN where
there is none), and `f_valid` is the fraction of each target cell covered by valid
source data. Threaded over rows.

With `nsub`, the source may be on any projection, and the mean is estimated from
`nsub` x `nsub` samples per target cell instead (see `remap_sampled.jl`).
"""
remap(tgt::ProjGrid, src::ProjGrid, F::AbstractMatrix; nsub::Union{Nothing,Integer}=nothing) =
    nsub === nothing ? remap(AlignedMap(tgt, src), F) : _remap_sampled(tgt, src, F, nsub)

function remap(m::AlignedMap, F::AbstractMatrix)
    size(F) == size(m.src) ||
        throw(DimensionMismatch("field size $(size(F)) does not match grid $(m.src.name) $(size(m.src))"))
    nxt, nyt = size(m.tgt)
    wx, wy = m.wx, m.wy

    # Source rows needed by any target row
    used = filter(!isempty, wy.ranges)
    ks = isempty(used) ? (1:0) : (minimum(first, used):maximum(last, used))
    koff = first(ks) - 1

    # Pass 1: along x, for each needed source row
    num = zeros(Float64, nxt, length(ks))
    den = zeros(Float64, nxt, length(ks))
    Threads.@threads for kk in eachindex(ks)
        k = ks[kk]
        @inbounds for i in 1:nxt
            a = 0.0
            b = 0.0
            for (n, l) in enumerate(wx.ranges[i])
                v = F[l, k]
                if _isvalid(v)
                    w = wx.weights[i][n]
                    a += w * v
                    b += w
                end
            end
            num[i, kk] = a
            den[i, kk] = b
        end
    end

    # Pass 2: along y, for each target row
    T = _outtype(F)
    Ft = Matrix{T}(undef, nxt, nyt)
    fv = Matrix{Float32}(undef, nxt, nyt)
    Threads.@threads for j in 1:nyt
        a = zeros(Float64, nxt)
        b = zeros(Float64, nxt)
        @inbounds for (n, k) in enumerate(wy.ranges[j])
            w = wy.weights[j][n]
            for i in 1:nxt
                a[i] += w * num[i, k-koff]
                b[i] += w * den[i, k-koff]
            end
        end
        @inbounds for i in 1:nxt
            Ft[i, j] = b[i] > 0 ? T(a[i] / b[i]) : T(NaN)
            fv[i, j] = Float32(b[i])
        end
    end
    return Ft, fv
end

"""
    remap_fractions(tgt::ProjGrid, src::ProjGrid, M, classes) -> (fracs, f_valid)

Area fraction of each class of the categorical field `M` within each target cell,
relative to its valid part. `fracs[c]` is the fraction of class `c`; `missing`
entries of `M` are excluded. With `nsub`, the source may be on any projection and
the fractions are estimated by supersampling.
"""
function remap_fractions(tgt::ProjGrid, src::ProjGrid, M::AbstractMatrix, classes;
                         nsub::Union{Nothing,Integer}=nothing)
    nsub === nothing || return _remap_fractions_sampled(tgt, src, M, classes, nsub)
    m = AlignedMap(tgt, src)
    fracs = Dict{eltype(classes),Matrix{Float32}}()
    fv = zeros(Float32, size(tgt))
    ind = Matrix{Float32}(undef, size(M))
    for c in classes
        Threads.@threads for k in axes(M, 2)
            @inbounds for l in axes(M, 1)
                v = M[l, k]
                ind[l, k] = ismissing(v) ? NaN32 : Float32(v == c)
            end
        end
        fracs[c], fv = remap(m, ind)
    end
    return fracs, fv
end

"""
    remap_dominant(tgt::ProjGrid, src::ProjGrid, M; parent=nothing) -> Mt
    remap_dominant(m::AlignedMap, M; parent=nothing) -> Mt

Class of the categorical field `M` (integers, size of `src`) covering the largest
area of each target cell, on the same projection (exact overlaps). Ties go to the
smaller class; target cells without source cells get 0. Threaded over rows.

With `parent = (Pt, Ps)`, a categorical field given on both grids (e.g. the coarser
level of a hierarchy of regions, already remapped), only the source cells with
`Ps == Pt[i, j]` count for target cell `(i, j)`, so the result stays nested in `Pt`.
Target cells where none of them overlaps get 0.
"""
remap_dominant(tgt::ProjGrid, src::ProjGrid, M::AbstractMatrix{<:Integer}; parent=nothing) =
    remap_dominant(AlignedMap(tgt, src), M; parent=parent)

function remap_dominant(m::AlignedMap, M::AbstractMatrix{T}; parent=nothing) where {T<:Integer}
    size(M) == size(m.src) ||
        throw(DimensionMismatch("field size $(size(M)) does not match grid $(m.src.name) $(size(m.src))"))
    if parent !== nothing
        Pt, Ps = parent
        (size(Pt) == size(m.tgt) && size(Ps) == size(m.src)) ||
            throw(DimensionMismatch("parent fields do not match the grids"))
    end
    nxt, nyt = size(m.tgt)
    wx, wy = m.wx, m.wy
    Mt = zeros(T, nxt, nyt)
    Threads.@threads for j in 1:nyt
        cls = T[]
        area = Float64[]
        @inbounds for i in 1:nxt
            empty!(cls)
            empty!(area)
            for (q, k) in enumerate(wy.ranges[j]), (p, l) in enumerate(wx.ranges[i])
                parent === nothing || Ps[l, k] == Pt[i, j] || continue
                w = wx.weights[i][p] * wy.weights[j][q]
                c = M[l, k]
                n = findfirst(==(c), cls)
                if n === nothing
                    push!(cls, c)
                    push!(area, w)
                else
                    area[n] += w
                end
            end
            isempty(cls) && continue
            best = 1
            for n in 2:length(cls)
                if area[n] > area[best] + 1e-12 || (abs(area[n] - area[best]) <= 1e-12 && cls[n] < cls[best])
                    best = n
                end
            end
            Mt[i, j] = cls[best]
        end
    end
    return Mt
end
