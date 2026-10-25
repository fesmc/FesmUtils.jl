"""
    FesmUtils.Namelist

Namelist (par) files as fesm-utils `nml.f90` reads them (branch `phase0-kryos`),
keeping the comments, and the configuration a run sees from them.

Files:

- `read_namelist(path)`: the groups of a file in file order (`Group` of `Key`),
  with every comment (before a group, after `&name`, before and after each key,
  after the last group).
- `write_namelist(path_or_io, groups; header="")`: writes the groups back
  (aligned keys and comments); reading the written file gives the same groups.
- `parse_value(raw)`, `format_value(v)`: the value of a key's text and back.
- `group_dict(groups)`: name => `Group`, the first copy winning (as `nml`).

Configuration:

- `Config(root; schema, library)`, `load_config(root)`: the defaults files
  (schema: which defaults group fills a group) and the library files (named
  groups not in the par file), from arguments or `.runme/info.json`.
- `resolve(cfg, base, overrides; follow)`: the groups of a base par file, the
  groups it names (`GroupRef` rules, `FORCING_REFS`) from the base, the library or
  the schema, and runme `-p` overrides (`parse_overrides`), every value with
  its source (`RGroup`, `Entry`, `value_of`).

Parsing rules (`nml.f90`):

- one key per line, `name = value`; the first `!` of a line starts a comment
  (also inside quotes, as `nml`); a comment after `&name` is allowed;
- a value is split into elements at blanks outside quotes ('...' or "..."),
  then quotes and a final comma are removed from each element:
  `"kg m-2 yr-1" "mm d-1"` is two elements, `-10, -5, 0` three;
- one element is a scalar, several a vector. An element that was quoted is a
  `String`; otherwise `Bool` (`T`, `True`, `.true.`, ...), `Int`, `Float64`
  (`d` exponents too) or `String`. A vector with a quoted element is a
  `Vector{String}`, mixed `Int`/`Float64` elements give `Vector{Float64}`.

Not supported: flat lines (`group.key = value`, an error) and continuation of
a value over several lines. Blanks are spaces and tabs.
"""
module Namelist

using JSON

export Key, Group, read_namelist, write_namelist, parse_value, format_value, group_dict
export Config, load_config, GroupRef, FORCING_REFS, Entry, RGroup, resolve, parse_overrides, value_of

"""
    Key(name, raw; comment="", pre="")

One key of a group: `raw` is the value as written (without the comment),
`comment` the inline comment (continuation lines joined by `\\n`), `pre` the
comment lines before the key inside the group.
"""
struct Key
    name::String
    raw::String
    comment::String
    pre::String
end
Key(name::AbstractString, raw::AbstractString; comment::AbstractString="", pre::AbstractString="") =
    Key(String(name), String(raw), String(comment), String(pre))

"""
    Group(name, keys; comment="", inline="", pre="", post="", after="")

A group `&name ... /`. Comments: `inline` after `&name` on its line, `comment`
the lines between `&name` and the first key, `post` the lines after the last
key, `pre` the lines before `&name` outside any group, `after` the lines after
`/` (only the last group of a file has them). Multi-line comments are joined by
`\\n`.
"""
struct Group
    name::String
    keys::Vector{Key}
    comment::String
    inline::String
    pre::String
    post::String
    after::String
end
Group(name::AbstractString, keys::AbstractVector{Key}=Key[]; comment::AbstractString="",
      inline::AbstractString="", pre::AbstractString="", post::AbstractString="",
      after::AbstractString="") =
    Group(String(name), collect(Key, keys), String(comment), String(inline), String(pre),
          String(post), String(after))

Base.:(==)(a::Key, b::Key) = a.name == b.name && a.raw == b.raw && a.comment == b.comment && a.pre == b.pre
Base.hash(k::Key, h::UInt) = hash((k.name, k.raw, k.comment, k.pre), h)
Base.:(==)(a::Group, b::Group) =
    a.name == b.name && a.keys == b.keys && a.comment == b.comment && a.inline == b.inline &&
    a.pre == b.pre && a.post == b.post && a.after == b.after
Base.hash(g::Group, h::UInt) = hash((g.name, g.keys, g.comment, g.inline, g.pre, g.post, g.after), h)

isblank(c::AbstractChar) = c == ' ' || c == '\t'

"Column (character index) of the first `!` of a line, or 0 (`nml` parse_line)."
function comment_column(line::AbstractString)
    for (i, c) in enumerate(line)
        c == '!' && return i
    end
    return 0
end

"The text before the first `!` and the comment after it, both stripped."
function split_comment(s::AbstractString)
    i = findfirst('!', s)
    i === nothing && return (String(strip(s)), "")
    return (String(strip(s[1:prevind(s, i)])), String(strip(s[nextind(s, i):end])))
end

"A flat line `group.key = value` (`nml` parse_flat_line, quote-aware)."
function is_flat(s::AbstractString)
    q = nothing
    dot = false
    for c in s
        if q !== nothing
            c == q && (q = nothing)
        elseif c == '"' || c == '\''
            q = c
        elseif c == '!'
            return false
        elseif c == '='
            return dot
        elseif c == '.'
            dot = true
        end
    end
    return false
end

"""
    read_namelist(path) -> Vector{Group}

The groups of a namelist file in file order, with their comments (see `Group`,
`Key`). A comment line inside a group at the column of the previous key's
inline `!` continues that comment. Blank lines are not kept.
"""
function read_namelist(path::AbstractString)
    groups = Group[]
    outside = String[]         # comment lines outside groups, for the next group
    name = nothing
    inline = ""
    gcomment = String[]
    keys = Key[]
    pre = String[]
    gpre = ""
    lastcol = 0                # column of the inline ! of the previous key
    for (ln, line) in enumerate(eachline(path))
        s = strip(line)
        if isempty(s)
            lastcol = 0
            continue
        end
        if name === nothing
            if startswith(s, "&")
                nm, inline = split_comment(s[2:end])
                isempty(nm) && error("$path:$ln: a group without a name")
                name = nm
                gpre = join(outside, "\n"); empty!(outside)
                gcomment = String[]; keys = Key[]; pre = String[]; lastcol = 0
            elseif startswith(s, "!")
                push!(outside, String(strip(s[2:end])))
            elseif is_flat(s)
                error("$path:$ln: flat lines (group.key = value) are not supported")
            end
            continue
        end
        if startswith(s, "/")
            push!(groups, Group(name, keys; comment=join(gcomment, "\n"), inline=inline,
                                pre=gpre, post=join(pre, "\n")))
            name = nothing
            continue
        end
        if startswith(s, "!")
            text = String(strip(s[2:end]))
            col = comment_column(line)
            if !isempty(keys) && lastcol > 0 && col == lastcol
                k = keys[end]
                keys[end] = Key(k.name, k.raw, k.comment * "\n" * text, k.pre)
            elseif isempty(keys)
                push!(gcomment, text)
            else
                push!(pre, text)
                lastcol = 0
            end
            continue
        end
        startswith(s, "&") && error("$path:$ln: group &$(strip(s[2:end])) inside &$name (missing /)")
        is_flat(s) && error("$path:$ln: flat lines (group.key = value) are not supported")
        occursin('=', s) || error("$path:$ln: no '=' on a parameter line of &$name")
        text, comment = split_comment(s)
        occursin('=', text) || error("$path:$ln: no '=' before the comment in &$name")
        k, v = split(text, '=', limit=2)
        push!(keys, Key(String(strip(k)), String(strip(v)), comment, join(pre, "\n")))
        pre = String[]
        lastcol = isempty(comment) ? 0 : comment_column(line)
    end
    name === nothing || error("$path: group &$name is not closed with /")
    if !isempty(outside)
        isempty(groups) && error("$path: comments but no group")
        g = groups[end]
        groups[end] = Group(g.name, g.keys; comment=g.comment, inline=g.inline, pre=g.pre,
                            post=g.post, after=join(outside, "\n"))
    end
    return groups
end

comment_lines(io::IO, text::AbstractString, indent::AbstractString) =
    isempty(text) || foreach(l -> println(io, indent, isempty(l) ? "!" : "! " * l), split(text, '\n'))

"""
    write_namelist(path_or_io, groups; header="")

Write groups as a namelist file: four-blank indent, `=` and inline comments
aligned per group, a blank line between groups. `header` (comment lines,
without `!`) goes before everything. `read_namelist` of the result gives back
`groups` (with `header` prepended to the first group's `pre`).
"""
function write_namelist(io::IO, groups::AbstractVector{Group}; header::AbstractString="")
    comment_lines(io, header, "")
    isempty(header) || isempty(groups) || println(io)
    for (n, g) in enumerate(groups)
        n > 1 && println(io)
        comment_lines(io, g.pre, "")
        println(io, "&", g.name, isempty(g.inline) ? "" : "   ! " * g.inline)
        comment_lines(io, g.comment, "    ")
        w = isempty(g.keys) ? 0 : maximum(k -> length(k.name), g.keys)
        heads = String[]
        for k in g.keys
            occursin('!', k.raw) && error("&$(g.name) $(k.name): a value with '!' cannot be written (nml reads a comment)")
            occursin('\n', k.raw) && error("&$(g.name) $(k.name): a value over several lines")
            push!(heads, "    " * rpad(k.name, w) * " = " * k.raw)
        end
        withc = [length(h) for (h, k) in zip(heads, g.keys) if !isempty(k.comment)]
        cw = isempty(withc) ? 0 : maximum(withc)
        for (h, k) in zip(heads, g.keys)
            comment_lines(io, k.pre, "    ")
            if isempty(k.comment)
                println(io, rstrip(h))
            else
                cl = split(k.comment, '\n')
                println(io, rpad(h, cw), "   ! ", cl[1])
                for l in cl[2:end]
                    println(io, " "^(cw + 3), isempty(l) ? "!" : "! " * l)
                end
            end
        end
        comment_lines(io, g.post, "    ")
        println(io, "/")
        if !isempty(g.after)
            println(io)
            comment_lines(io, g.after, "")
        end
    end
    return nothing
end
write_namelist(path::AbstractString, groups::AbstractVector{Group}; header::AbstractString="") =
    open(io -> write_namelist(io, groups; header=header), path, "w")

"Groups as a Dict name => Group (the first copy wins, as `nml`)."
function group_dict(groups::AbstractVector{Group})
    d = Dict{String,Group}()
    for g in groups
        haskey(d, g.name) || (d[g.name] = g)
    end
    return d
end

"Elements of a value: split at blanks outside quotes (`nml` string_to_vector), as (text, quoted)."
function elements(raw::AbstractString)
    out = Tuple{String,Bool}[]
    buf = IOBuffer()
    q = nothing
    quoted = false
    started = false
    function push_element!()
        t = String(take!(buf))
        t = replace(t, '"' => "", '\'' => "")      # remove_quotes_comma: every quote
        t = String(strip(t))
        endswith(t, ',') && (t = String(strip(t[1:prevind(t, lastindex(t))])))
        push!(out, (t, quoted))
        quoted = false
        started = false
    end
    for c in raw
        if q !== nothing
            print(buf, c)
            c == q && (q = nothing)
        elseif isblank(c)
            started && push_element!()
        else
            started = true
            print(buf, c)
            if c == '"' || c == '\''
                q = c
                quoted = true
            end
        end
    end
    started && push_element!()
    return out
end

const TRUE_WORDS = ("T", "True", "TRUE", "true", ".TRUE.", ".true.")
const FALSE_WORDS = ("F", "False", "FALSE", "false", ".FALSE.", ".false.")

function parse_scalar(t::AbstractString, quoted::Bool)
    quoted && return String(t)
    t in TRUE_WORDS && return true
    t in FALSE_WORDS && return false
    occursin(r"^[-+]?\d+$", t) && return parse(Int, t)
    if occursin(r"^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$", t)
        return parse(Float64, replace(t, r"[dD]" => "e"))
    end
    return String(t)
end

"""
    parse_value(raw)

The Julia value of a key's text (see the module docstring): a `String`,
`Bool`, `Int` or `Float64`, or a `Vector` of them. `""` and an empty text give
the empty string.
"""
function parse_value(raw::AbstractString)
    els = elements(raw)
    isempty(els) && return ""
    vals = [parse_scalar(t, q) for (t, q) in els]
    length(vals) == 1 && return vals[1]
    any(last, els) && return [t for (t, _) in els]
    all(v -> v isa Bool, vals) && return Bool[v for v in vals]
    all(v -> v isa Int, vals) && return Int[v for v in vals]
    all(v -> v isa Real && !(v isa Bool), vals) && return Float64[v for v in vals]
    return [t for (t, _) in els]
end

format_scalar(v::Bool) = v ? "True" : "False"
format_scalar(v::Integer) = string(v)
function format_scalar(v::AbstractFloat)
    isfinite(v) || error("format_value: $v cannot be written to a namelist")
    return string(Float64(v))
end
function format_scalar(v::AbstractString)
    (occursin('"', v) || occursin('\'', v)) && error("format_value: a string with a quote: $v")
    return "\"" * v * "\""
end

"""
    format_value(v)

The text of a value as a par file writes it: strings quoted, `True`/`False`,
floats in the shortest form that reads back to the same `Float64`, vectors
joined by blanks. `parse_value(format_value(v)) == v` for the values
`parse_value` returns (a vector of one element reads back as a scalar).
"""
format_value(v) = format_scalar(v)
format_value(v::AbstractVector) = join(format_scalar.(v), " ")

include("namelist_resolve.jl")

end # module
