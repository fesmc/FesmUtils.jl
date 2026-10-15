# Exact conservative remapping onto a projected grid from a regular lon-lat grid.
#
# Each source cell is a quadrilateral on the target projection, from its four
# projected corners (its edges taken as straight: for cells of a few arc seconds the
# area is then exact to about 1e-9). It is clipped against each target cell it
# overlaps, and the area of the clipped polygon is its weight in that cell, as a
# fraction of the target cell area. Unlike sampling (`SampledMap`), every part of a
# source cell counts in proportion to its area, so fields that vary at the scale of
# the source cells are not aliased, also on target cells smaller than them.

"""
    ConservativeMap(tgt::ProjGrid, src::LonLatGrid)
    remap(m::ConservativeMap, F) -> (Ft, f_valid)

Exact conservative remapping from the lon-lat grid `src` onto the projected grid
`tgt`: for each target cell, the source cells that overlap it and the fraction of its
area that each covers. Build it once and reuse it for many fields on the same grids.
Only the source cells over the target are projected; a global source may wrap in
longitude across the target. Its size grows with the number of overlapping pairs, so
build it for sources not much finer than needed. Threaded, and reproducible (the
weights of each target cell are in source order).
"""
struct ConservativeMap
    tgt::ProjGrid
    src::LonLatGrid
    ptr::Vector{Int}        # weights of target cell k (linear index): ptr[k]:ptr[k+1]-1
    idx::Vector{Int}        # linear index of the source cell
    w::Vector{Float64}      # fraction of the target cell area that it covers
end

# Columns (unwrapped: index k is longitude edge lon[1] - dlon/2 + (k-1) dlon, column
# mod1(k, n)) and rows of the source cells that may overlap the projected grid `g`,
# from its outline, with a margin of 2 cells
function _source_window(g::ProjGrid, src::LonLatGrid; npts::Integer=200)
    nlon, nlat = length(src.lon), length(src.lat)
    latmin, latmax = lat_bounds(g; margin=0.0)
    j0 = max(1, searchsortedlast(src.latedges, latmin) - 2)
    j1 = min(nlat, searchsortedfirst(src.latedges, latmax) + 2)
    j0 <= j1 || return 1:0, 1:0
    if latmin <= -90 || latmax >= 90            # a pole inside: all longitudes
        return 1:nlon, j0:j1
    end
    dx, dy = spacing(g)
    x0, x1 = g.xc[1] - dx / 2, g.xc[end] + dx / 2
    y0, y1 = g.yc[1] - dy / 2, g.yc[end] + dy / 2
    s = range(0, 1; length=npts)
    pts = vcat([(x0 + t * (x1 - x0), y) for t in s, y in (y0, y1)][:],
               [(x, y0 + t * (y1 - y0)) for t in s, x in (x0, x1)][:])
    trans = _inverse_transform(g)
    lon = first.(map(trans, pts))
    ref = lon[1]
    lon = [ref + mod(l - ref + 180, 360) - 180 for l in lon]   # unwrapped around the outline
    lonmin, lonmax = extrema(lon)
    e1 = src.lon[1] - src.dlon / 2
    if src.isglobal
        lonmin = e1 + mod(lonmin - e1, 360)
        lonmax = lonmin + (maximum(lon) - minimum(lon))
        i0 = floor(Int, (lonmin - e1) / src.dlon) + 1 - 2
        i1 = floor(Int, (lonmax - e1) / src.dlon) + 1 + 2
        return (i1 - i0 + 1 >= nlon ? (1:nlon) : (i0:i1)), j0:j1
    end
    i0 = max(1, floor(Int, (lonmin - e1) / src.dlon) + 1 - 2)
    i1 = min(nlon, floor(Int, (lonmax - e1) / src.dlon) + 1 + 2)
    return i0:i1, j0:j1
end

# Area of the polygon (px, py) (shoelace; positive if counter-clockwise)
function _polygon_area(px, py, n)
    a = 0.0
    @inbounds for k in 1:n
        l = k == n ? 1 : k + 1
        a += px[k] * py[l] - px[l] * py[k]
    end
    return a / 2
end

# Clip the polygon (px, py, n points) to the half-plane s * (c - bound) >= 0 of
# coordinate c (x if onx, else y), into (qx, qy); returns the number of points
@inline function _clip_halfplane!(qx, qy, px, py, n, onx::Bool, bound::Float64, s::Float64)
    m = 0
    n == 0 && return 0
    @inbounds for k in 1:n
        l = k == n ? 1 : k + 1
        ak = s * ((onx ? px[k] : py[k]) - bound)
        al = s * ((onx ? px[l] : py[l]) - bound)
        if ak >= 0
            m += 1; qx[m] = px[k]; qy[m] = py[k]
        end
        if (ak >= 0) != (al >= 0)
            t = ak / (ak - al)
            m += 1; qx[m] = px[k] + t * (px[l] - px[k]); qy[m] = py[k] + t * (py[l] - py[k])
        end
    end
    return m
end

# Area of the quadrilateral (px, py) within the rectangle [xa, xb] x [ya, yb]
# (Sutherland-Hodgman), with work buffers a and b (8 points each)
function _clipped_area(px, py, xa, xb, ya, yb, ax, ay, bx, by)
    n = _clip_halfplane!(ax, ay, px, py, 4, true, xa, 1.0)
    n = _clip_halfplane!(bx, by, ax, ay, n, true, xb, -1.0)
    n = _clip_halfplane!(ax, ay, bx, by, n, false, ya, 1.0)
    n = _clip_halfplane!(bx, by, ax, ay, n, false, yb, -1.0)
    return n < 3 ? 0.0 : abs(_polygon_area(bx, by, n))
end

function ConservativeMap(tgt::ProjGrid, src::LonLatGrid)
    nx, ny = size(tgt)
    nlon = length(src.lon)
    dx, dy = spacing(tgt)
    xe0, ye0 = tgt.xc[1] - dx / 2, tgt.yc[1] - dy / 2
    cols, rows = _source_window(tgt, src)
    nc, nr = length(cols), length(rows)

    # Projected corners of the source cells in the window (threaded over rows)
    e1 = src.lon[1] - src.dlon / 2
    lone = [e1 + (k - 1) * src.dlon for k in first(cols):last(cols)+1]
    late = src.latedges[first(rows):last(rows)+1]
    cx = Matrix{Float64}(undef, nc + 1, nr + 1)
    cy = Matrix{Float64}(undef, nc + 1, nr + 1)
    if nc > 0 && nr > 0
        tasks = map(_chunks(nr + 1)) do js
            Threads.@spawn begin
                trans = _forward_transform(tgt)
                for j in js, i in 1:nc+1
                    cx[i, j], cy[i, j] = trans((lone[i], late[j]))
                end
            end
        end
        foreach(wait, tasks)
    end

    # Overlaps of each source cell with the target cells (threaded over source rows),
    # as (target, source, area) per chunk
    parts = map(collect(_chunks(nr))) do js
        Threads.@spawn begin
            out = Tuple{Int,Int,Float64}[]
            px, py, qx, qy = zeros(4), zeros(4), zeros(4), zeros(4)
            ax, ay, bx, by = zeros(8), zeros(8), zeros(8), zeros(8)
            for jj in js, ii in 1:nc
                px[1], py[1] = cx[ii, jj], cy[ii, jj]
                px[2], py[2] = cx[ii+1, jj], cy[ii+1, jj]
                px[3], py[3] = cx[ii+1, jj+1], cy[ii+1, jj+1]
                px[4], py[4] = cx[ii, jj+1], cy[ii, jj+1]
                all(isfinite, px) && all(isfinite, py) || continue
                xmin, xmax = extrema(px)
                ymin, ymax = extrema(py)
                i0 = max(1, floor(Int, (xmin - xe0) / dx) + 1)
                i1 = min(nx, floor(Int, (xmax - xe0) / dx) + 1)
                j0 = max(1, floor(Int, (ymin - ye0) / dy) + 1)
                j1 = min(ny, floor(Int, (ymax - ye0) / dy) + 1)
                (i0 <= i1 && j0 <= j1) || continue
                s = mod1(cols[ii], nlon) + (rows[jj] - 1) * nlon
                for j in j0:j1, i in i0:i1
                    # In coordinates relative to the target cell (precise areas)
                    xa = xe0 + (i - 1) * dx
                    ya = ye0 + (j - 1) * dy
                    for k in 1:4
                        qx[k] = px[k] - xa
                        qy[k] = py[k] - ya
                    end
                    a = _clipped_area(qx, qy, 0.0, dx, 0.0, dy, ax, ay, bx, by)
                    a > 1e-12 * dx * dy && push!(out, (i + (j - 1) * nx, s, a / (dx * dy)))
                end
            end
            out
        end
    end
    pairs = reduce(vcat, fetch.(parts); init=Tuple{Int,Int,Float64}[])
    sort!(pairs; by=p -> (p[1], p[2]))

    # Compressed by target cell; a source cell met twice (a global source wrapping
    # across the window) is summed
    ptr = Vector{Int}(undef, nx * ny + 1)
    idx = Int[]
    w = Float64[]
    sizehint!(idx, length(pairs)); sizehint!(w, length(pairs))
    p = 1
    ptr[1] = 1
    for k in 1:nx*ny
        while p <= length(pairs) && pairs[p][1] == k
            s, a = pairs[p][2], 0.0
            while p <= length(pairs) && pairs[p][1] == k && pairs[p][2] == s
                a += pairs[p][3]
                p += 1
            end
            push!(idx, s)
            push!(w, a)
        end
        ptr[k+1] = length(idx) + 1
    end
    return ConservativeMap(tgt, src, ptr, idx, w)
end

function remap(m::ConservativeMap, F::AbstractMatrix)
    size(F) == size(m.src) ||
        throw(DimensionMismatch("field size $(size(F)) does not match source grid $(size(m.src))"))
    nx, ny = size(m.tgt)
    T = _outtype(F)
    Ft = Matrix{T}(undef, nx, ny)
    fv = Matrix{Float32}(undef, nx, ny)
    Threads.@threads for j in 1:ny
        @inbounds for i in 1:nx
            k = i + (j - 1) * nx
            num, a = 0.0, 0.0
            for p in m.ptr[k]:m.ptr[k+1]-1
                v = F[m.idx[p]]
                _isvalid(v) || continue
                num += m.w[p] * v
                a += m.w[p]
            end
            Ft[i, j] = a > 0 ? T(num / a) : T(NaN)
            fv[i, j] = Float32(a)
        end
    end
    return Ft, fv
end
