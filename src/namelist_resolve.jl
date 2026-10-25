# The configuration a run sees: a base par file, the groups it names from the
# library files, the defaults (schema) groups and the runme -p overrides, with
# the source of every value. Part of module Namelist.

"""
    Entry

One resolved key: `raw` text, parsed `value`, `source` (`base`, `library`,
`default` or `override`), `comment` and `pre` (from the defaults group when
the par file has none).
"""
struct Entry
    raw::String
    value::Any
    source::String
    comment::String
    pre::String
end

"""
    RGroup

One resolved group: `name`, `source` (`base`, `library` or `new`: only
defaults), `file` it was read from (`""` for new), `template` (defaults file,
defaults group) or `nothing`, the key `order` and the `entries`
(name => `Entry`). `Group(rg)` gives a writable group.
"""
mutable struct RGroup
    name::String
    source::String
    file::String
    template::Union{Nothing,Tuple{String,String}}
    order::Vector{String}
    entries::Dict{String,Entry}
    comment::String
end

"`value_of(rg, key, default=nothing)`: the value of a key of a resolved group."
value_of(rg::RGroup, key::AbstractString, default=nothing) =
    haskey(rg.entries, key) ? rg.entries[key].value : default

Group(rg::RGroup) = Group(rg.name, [Key(k, rg.entries[k].raw, rg.entries[k].comment, rg.entries[k].pre)
                                    for k in rg.order]; comment=rg.comment)

"""
    GroupRef(from, key; prefix="", when=nothing, template=nothing)

A group named by a key: in every resolved group whose name matches `from`
(a `Regex`), each name listed in `key` (`""` entries skipped) names the group
`prefix * name`, which is resolved too. `{k}` in `key` is replaced by the value
of key `k` of the same group (`"nml_{model}"`). `when = "k" => v`: only when key
`k` is `v` (an absent `k` counts as `v`). `template = (file, group)`: the
defaults group of the named group (as `nml_read(..., defaults_group)`).
"""
struct GroupRef
    from::Regex
    key::String
    prefix::String
    when::Union{Nothing,Pair{String,<:Any}}
    template::Union{Nothing,Tuple{String,String}}
end
GroupRef(from::Regex, key::AbstractString; prefix::AbstractString="", when=nothing, template=nothing) =
    GroupRef(from, String(key), String(prefix), when, template)

"""
The named groups of the fesm-utils forcing (tsgen coefficients, climate terms,
snapshots): `[forcing<sfx>]` `coefs`, `reference`, `terms`, and `[term_<name>]`
`snap0`, `snap1` of snapshot terms.
"""
const FORCING_REFS = [
    GroupRef(r"^forcing(_\w+)?$", "coefs"; prefix="coef_"),
    GroupRef(r"^forcing(_\w+)?$", "reference"; prefix="snapshot_"),
    GroupRef(r"^forcing(_\w+)?$", "terms"; prefix="term_"),
    GroupRef(r"^term_\w+$", "snap0"; prefix="snapshot_", when="pattern" => "snapshots"),
    GroupRef(r"^term_\w+$", "snap1"; prefix="snapshot_", when="pattern" => "snapshots"),
]

"""
    Config(root; schema=[], library=String[])

Where the groups of a run come from, paths relative to `root`: `schema`, a
list of `name => (defaults_file, defaults_group)` (`name` a `String` or a
`Regex`; the first match wins) giving the defaults of a group, and `library`,
the files searched in order for a named group that is not in the par file.
"""
struct Config
    root::String
    schema::Vector{Pair{Any,Tuple{String,String}}}
    library::Vector{String}
    files::Dict{String,Vector{Group}}
end
Config(root::AbstractString; schema=Pair{Any,Tuple{String,String}}[], library=String[]) =
    Config(String(root), Pair{Any,Tuple{String,String}}[p for p in schema], String.(library),
           Dict{String,Vector{Group}}())

"Absolute path of a file of the configuration."
cfg_path(cfg::Config, f::AbstractString) = isabspath(f) ? String(f) : joinpath(cfg.root, f)

"The groups of a file (read once)."
cfg_groups(cfg::Config, f::AbstractString) = get!(() -> read_namelist(cfg_path(cfg, f)), cfg.files, f)

function cfg_find(cfg::Config, file::AbstractString, group::AbstractString)
    for g in cfg_groups(cfg, file)
        g.name == group && return g
    end
    return nothing
end

"The defaults (file, group) of a group name, or `nothing`."
function template(cfg::Config, name::AbstractString)
    for (pat, t) in cfg.schema
        (pat isa Regex ? occursin(pat, name) : pat == name) && return t
    end
    return nothing
end

"""
The groups named like the template groups of fesm-utils' forcing, as
`Regex => defaults group` (`load_config`).
"""
const FORCING_TEMPLATES = [r"^coef_\w+$" => "coef", r"^term_\w+$" => "term", r"^snapshot_\w+$" => "snapshot"]

"""
    load_config(root; templates=FORCING_TEMPLATES)

The `Config` of a runme repository from `.runme/info.json`: every group of the
`par_defaults` files is the defaults of the group of the same name (the first
file wins), `templates` (`Regex => group`) name the defaults of named groups
(found in the first defaults file that has the group), and `par_library` are
the library files.
"""
function load_config(root::AbstractString; templates=FORCING_TEMPLATES)
    info = JSON.parsefile(joinpath(root, ".runme", "info.json"))
    dfiles = String.(get(info, "par_defaults", String[]))
    lfiles = String.(get(info, "par_library", String[]))
    cfg = Config(root; library=lfiles)
    seen = Set{String}()
    for f in dfiles, g in cfg_groups(cfg, f)
        g.name in seen && continue
        push!(seen, g.name)
        push!(cfg.schema, g.name => (f, g.name))
    end
    for (re, t) in templates
        f = findfirst(f -> cfg_find(cfg, f, t) !== nothing, dfiles)
        f === nothing || push!(cfg.schema, re => (dfiles[f], t))
    end
    return cfg
end

"""
    parse_overrides(args) -> Vector{Tuple{String,Vector{String},Bool}}

The `-p` arguments of runme as `(group.key, tokens, vector)`: `g.k=v` a
scalar, `g.k==a,b` a vector. Ensemble forms (a comma list after a single `=`)
are rejected.
"""
function parse_overrides(args::AbstractVector{<:AbstractString})
    out = Tuple{String,Vector{String},Bool}[]
    for a in args
        if occursin("==", a)
            k, v = split(a, "==", limit=2)
            push!(out, (String(k), String.(split(v, ",")), true))
        else
            occursin("=", a) || error("override '$a': expected group.key=value")
            k, v = split(a, "=", limit=2)
            occursin(",", v) && error("override '$a': an ensemble list; use == for a vector")
            push!(out, (String(k), [String(v)], false))
        end
    end
    return out
end

"Coerce the text of an override to the type of the value it replaces."
function coerce(tok::AbstractString, ref)
    t = String(strip(tok, ['"', '\'']))
    if ref isa Bool
        t in TRUE_WORDS && return true
        t in FALSE_WORDS && return false
        error("override value '$tok': expected a logical")
    end
    if ref isa Integer
        x = tryparse(Int, t)
        x === nothing && error("override value '$tok': expected an integer")
        return x
    end
    if ref isa Real
        x = tryparse(Float64, replace(t, r"[dD]" => "e"))
        x === nothing && error("override value '$tok': expected a number")
        return x
    end
    return t
end

function override_value(toks, vector::Bool, ref)
    refel = ref isa AbstractVector ? (isempty(ref) ? "" : ref[1]) : ref
    vals = [coerce(t, refel) for t in toks]
    (vector || ref isa AbstractVector) && length(vals) > 1 && return vals
    return vals[1]
end

"Names listed in a key (`\"\"` entries dropped)."
function names_of(rg::RGroup, key::AbstractString)
    v = value_of(rg, key, "")
    v isa AbstractVector && return [string(x) for x in v if !isempty(string(x))]
    s = string(v)
    return isempty(s) ? String[] : [s]
end

"`{k}` in a key name replaced by the value of key k of the group (nothing if absent)."
function ref_key(rg::RGroup, key::AbstractString)
    out = key
    for m in eachmatch(r"\{(\w+)\}", key)
        v = value_of(rg, m.captures[1])
        (v === nothing || v isa AbstractVector) && return nothing
        out = replace(out, m.match => string(v))
    end
    return out
end

"""
    resolve(cfg, base, overrides=[]; follow=FORCING_REFS) -> (groups, order)

Resolve the par file `base` (a path): its groups, then the groups named by the
`follow` rules (`GroupRef`) from the base, else the library files of `cfg`, else
new from their defaults group, then the overrides (`parse_overrides`), then
the named groups again. A group with a defaults group (`cfg` schema, or the
`template` of the `GroupRef` naming it) gets every key of it, the par file's
values over the defaults. Returns `Dict` name => `RGroup` and the group names
in order (the base's, then the named ones as found). Errors name the group or
key that cannot be resolved.
"""
function resolve(cfg::Config, base::AbstractString, overrides=Tuple{String,Vector{String},Bool}[];
                 follow::AbstractVector{GroupRef}=FORCING_REFS)
    base_groups = read_namelist(base)
    based = group_dict(base_groups)
    out = Dict{String,RGroup}()
    order = String[]

    function fill_defaults!(rg::RGroup, t::Tuple{String,String})
        dg = cfg_find(cfg, t[1], t[2])
        dg === nothing && error("defaults group &$(t[2]) of &$(rg.name) is not in $(t[1])")
        rg.template = t
        neworder = String[]
        for k in dg.keys
            push!(neworder, k.name)
            if haskey(rg.entries, k.name)
                e = rg.entries[k.name]
                rg.entries[k.name] = Entry(e.raw, e.value, e.source,
                                           isempty(e.comment) ? k.comment : e.comment,
                                           isempty(e.pre) ? k.pre : e.pre)
            else
                rg.entries[k.name] = Entry(k.raw, parse_value(k.raw), "default", k.comment, k.pre)
            end
        end
        append!(neworder, filter(k -> !(k in neworder), rg.order))
        rg.order = neworder
        isempty(rg.comment) && (rg.comment = dg.comment)
        return rg
    end

    function group!(name::AbstractString, t=nothing)
        if haskey(out, name)
            rg = out[name]
            t !== nothing && rg.template === nothing && fill_defaults!(rg, t)
            return rg
        end
        t === nothing && (t = template(cfg, name))
        if haskey(based, name)
            g = based[name]; src = "base"; file = String(base)
        else
            g = nothing; src = "new"; file = ""
            for f in cfg.library
                g = cfg_find(cfg, f, name)
                g === nothing || (src = "library"; file = cfg_path(cfg, f); break)
            end
            if g === nothing
                t === nothing && error("group &$name: not in the par file, the library ($(join(cfg.library, ", "))) or the defaults")
                g = Group(name)
            end
        end
        entries = Dict{String,Entry}()
        ord = String[]
        for k in g.keys
            haskey(entries, k.name) && continue     # nml: the first line wins
            entries[k.name] = Entry(k.raw, parse_value(k.raw), src, k.comment, k.pre)
            push!(ord, k.name)
        end
        rg = RGroup(String(name), src, file, nothing, ord, entries, g.comment)
        t === nothing || fill_defaults!(rg, t)
        out[name] = rg
        push!(order, name)
        return rg
    end

    function follow!()
        changed = true
        while changed
            changed = false
            for gname in copy(order), r in follow
                occursin(r.from, gname) || continue
                rg = out[gname]
                if r.when !== nothing
                    v = value_of(rg, r.when.first, r.when.second)
                    v == r.when.second || continue
                end
                key = ref_key(rg, r.key)
                key === nothing && continue
                for n in names_of(rg, key)
                    name = r.prefix * n
                    before = (haskey(out, name), haskey(out, name) && out[name].template !== nothing)
                    group!(name, r.template)
                    after = (true, out[name].template !== nothing)
                    before == after || (changed = true)
                end
            end
        end
    end

    for g in base_groups
        group!(g.name)
    end
    follow!()
    for (key, toks, vector) in overrides
        occursin(".", key) || error("override '$key': expected group.key")
        gname, kname = split(key, ".", limit=2)
        rg = group!(gname)
        haskey(rg.entries, kname) || error("override '$key': key not in &$gname nor its defaults")
        e = rg.entries[kname]
        v = override_value(toks, vector, e.value)
        rg.entries[kname] = Entry(format_value(v), v, "override", e.comment, e.pre)
    end
    follow!()
    return out, order
end
