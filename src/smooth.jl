# Gaussian smoothing of fields with missing values.
#
# The field is convolved with a Gaussian kernel (truncated at 3 standard deviations)
# in two separable passes, along x (or longitude) and then along y (or latitude).
# Missing values (`missing` or NaN) are excluded by normalised convolution: the
# convolution of the valid values is divided by that of the valid-cell indicator, so
# each valid cell gets a weighted mean over its valid neighbours, and missing cells
# stay missing.

# Gaussian weights at offsets -r:r cells, for a standard deviation of `s` cells.
function _gauss_weights(s::Real)
    s > 0 || return [1.0]
    r = ceil(Int, 3s)
    return [exp(-0.5 * (k / s)^2) for k in -r:r]
end

# One pass along the first dimension of (num, den), in place: row j is convolved with
# weights(j) (odd length, centred), periodically if `periodic`.
function _smooth_pass_x!(num::Matrix{Float64}, den::Matrix{Float64}, weights, periodic::Bool)
    nx, ny = size(num)
    Threads.@threads for j in 1:ny
        w = weights(j)
        r = length(w) ÷ 2
        length(w) == 1 && continue
        a = num[:, j]
        b = den[:, j]
        @inbounds for i in 1:nx
            sa = 0.0
            sb = 0.0
            for k in -r:r
                l = i + k
                if periodic
                    l = mod1(l, nx)
                elseif !(1 <= l <= nx)
                    continue
                end
                wk = w[k+r+1]
                sa += wk * a[l]
                sb += wk * b[l]
            end
            num[i, j] = sa
            den[i, j] = sb
        end
    end
end

# One pass along the second dimension, with weights w (odd length, centred) times the
# row weights `rowweight` (e.g. the cell areas of a lon-lat grid).
function _smooth_pass_y(num::Matrix{Float64}, den::Matrix{Float64}, w::Vector{Float64}, rowweight)
    nx, ny = size(num)
    r = length(w) ÷ 2
    outn = zeros(Float64, nx, ny)
    outd = zeros(Float64, nx, ny)
    Threads.@threads for j in 1:ny
        @inbounds for k in max(1, j - r):min(ny, j + r)
            wk = w[k-j+r+1] * rowweight[k]
            for i in 1:nx
                outn[i, j] += wk * num[i, k]
                outd[i, j] += wk * den[i, k]
            end
        end
    end
    return outn, outd
end

function _smooth(F::AbstractMatrix, xweights, yweights::Vector{Float64}, rowweight, periodic::Bool)
    valid = map(_isvalid, F)
    num = map((v, ok) -> ok ? Float64(v) : 0.0, F, valid)
    den = Float64.(valid)
    _smooth_pass_x!(num, den, xweights, periodic)
    num, den = _smooth_pass_y(num, den, yweights, rowweight)
    T = _outtype(F)
    return map((a, b, ok) -> ok && b > 0 ? T(a / b) : T(NaN), num, den, valid)
end

"""
    smooth(g::ProjGrid, F, sigma) -> Fs
    smooth(g::LonLatGrid, F, sigma) -> Fs

Gaussian smoothing of the field `F` (size of `g`) with standard deviation `sigma`
(km). Missing values (`missing` or NaN) are excluded: each valid cell gets the
weighted mean of its valid neighbours, and missing cells stay missing. On a lon-lat
grid, the kernel is `sigma` km wide in both directions at every latitude (wider in
longitude towards the poles, up to the whole circle), periodic in longitude on a
global grid, and weighted by cell area. Separable, threaded.
"""
function smooth(g::ProjGrid, F::AbstractMatrix, sigma::Real)
    size(F) == size(g) || throw(DimensionMismatch("field size $(size(F)) does not match grid $(size(g))"))
    sigma >= 0 || throw(ArgumentError("sigma must be >= 0"))
    dx, dy = spacing(g)
    wx = _gauss_weights(sigma / dx)
    return _smooth(F, _ -> wx, _gauss_weights(sigma / dy), ones(size(g)[2]), false)
end

function smooth(g::LonLatGrid, F::AbstractMatrix, sigma::Real)
    size(F) == size(g) || throw(DimensionMismatch("field size $(size(F)) does not match grid $(size(g))"))
    sigma >= 0 || throw(ArgumentError("sigma must be >= 0"))
    km = 111.195                                   # km per degree on the mean sphere
    nx = length(g.lon)
    h = (nx - 1) ÷ 2                               # widest kernel: each longitude once
    wx = map(g.lat) do lat
        w = _gauss_weights(sigma / (g.dlon * km * max(cosd(lat), 1e-6)))
        r = length(w) ÷ 2
        r <= h ? w : w[r+1-h:r+1+h]
    end
    wy = _gauss_weights(sigma / (g.dlat * km))
    return _smooth(F, j -> wx[j], wy, cosd.(g.lat), g.isglobal)
end
