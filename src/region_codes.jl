# Hierarchical region codes.
#
# A region is named by its path in a tree of regions, e.g. "1.3.1" for the first
# subregion of region 3 of region 1. Its code is the integer with two decimal
# digits per level below the first: "1" => 1, "1.3" => 103, "1.3.1" => 10301. Each
# level has at most 99 regions, and the level of a code follows from its size. A
# field of codes for level n holds the code of the deepest region at or above level
# n that covers each cell, so cells without a subregion keep the code of their
# parent region.

const _REGION_BASE = 100

"""
    region_code(path) -> Int

Code of the region with path `path` ("1.3.1" => 10301).
"""
function region_code(path::AbstractString)
    parts = parse.(Int, split(path, '.'))
    all(p -> 1 <= p < _REGION_BASE, parts) ||
        throw(ArgumentError("region path $path: each part must be in 1..$(_REGION_BASE - 1)"))
    return foldl((c, p) -> c * _REGION_BASE + p, parts)
end

"""
    region_level(code) -> Int

Level of a region code (1 for 1..99, 2 for 101..9999, ...).
"""
function region_level(code::Integer)
    code >= 1 || throw(ArgumentError("invalid region code $code"))
    n = 1
    while code >= _REGION_BASE
        code ÷= _REGION_BASE
        n += 1
    end
    return n
end

"""
    region_path(code) -> String

Path of a region code (10301 => "1.3.1").
"""
function region_path(code::Integer)
    parts = Int[]
    c = Int(code)
    for _ in 1:region_level(code)
        pushfirst!(parts, c % _REGION_BASE)
        c ÷= _REGION_BASE
    end
    any(iszero, parts) && throw(ArgumentError("invalid region code $code"))
    return join(parts, '.')
end

"""
    region_ancestor(code, level) -> Int

Code of the region at `level` that contains region `code` (`code` itself at its own
level).
"""
function region_ancestor(code::Integer, level::Integer)
    n = region_level(code)
    1 <= level <= n || throw(ArgumentError("region $code has no ancestor at level $level"))
    return code ÷ _REGION_BASE^(n - level)
end

"""
    in_region(codes, code) -> BitArray

Cells of the field of region codes `codes` that lie in region `code` or one of its
subregions. Cells with code 0 are in no region.
"""
function in_region(codes::AbstractArray{<:Integer}, code::Integer)
    level = region_level(code)
    return map(codes) do c
        c > 0 && region_level(c) >= level && region_ancestor(c, level) == code
    end
end

"""
    region_flag_attrib(codes, names) -> Vector{Pair{String,Any}}

CF attributes `flag_values` and `flag_meanings` for a field of region codes, with
`names` the region names (spaces become underscores).
"""
function region_flag_attrib(codes::AbstractVector{<:Integer}, names::AbstractVector{<:AbstractString})
    length(codes) == length(names) || throw(ArgumentError("codes and names differ in length"))
    return Pair{String,Any}["flag_values" => Int32.(codes),
                            "flag_meanings" => join(replace.(names, r"\s+" => "_"), " ")]
end
