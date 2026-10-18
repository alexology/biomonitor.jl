# Minimal JSON parser for the responses of the taxonomic web services.
# Objects become Dict{String,Any}, arrays Vector{Any}, null becomes nothing.

struct _JSONReader
    s::String
    i::Base.RefValue{Int}
end

parse_json(s::AbstractString) = (r = _JSONReader(String(s), Ref(1)); v = _jvalue(r); _jskip(r);
                                 r.i[] <= ncodeunits(r.s) && throw(ArgumentError("trailing characters in JSON")); v)

function _jskip(r)
    while r.i[] <= ncodeunits(r.s) && r.s[r.i[]] in (' ', '\t', '\n', '\r')
        r.i[] = nextind(r.s, r.i[])
    end
end

_jpeek(r) = (_jskip(r); r.i[] <= ncodeunits(r.s) ? r.s[r.i[]] : throw(ArgumentError("unexpected end of JSON")))
_jnext!(r) = (c = r.s[r.i[]]; r.i[] = nextind(r.s, r.i[]); c)

function _jexpect(r, word)
    for c in word
        _jnext!(r) == c || throw(ArgumentError("invalid JSON literal"))
    end
end

function _jvalue(r)
    c = _jpeek(r)
    if c == '{'
        _jnext!(r)
        d = Dict{String,Any}()
        _jpeek(r) == '}' && (_jnext!(r); return d)
        while true
            k = _jstring(r)
            _jpeek(r) == ':' || throw(ArgumentError("expected : in JSON object"))
            _jnext!(r)
            d[k] = _jvalue(r)
            c2 = _jpeek(r)
            _jnext!(r)
            c2 == '}' && return d
            c2 == ',' || throw(ArgumentError("expected , or } in JSON object"))
        end
    elseif c == '['
        _jnext!(r)
        v = Any[]
        _jpeek(r) == ']' && (_jnext!(r); return v)
        while true
            push!(v, _jvalue(r))
            c2 = _jpeek(r)
            _jnext!(r)
            c2 == ']' && return v
            c2 == ',' || throw(ArgumentError("expected , or ] in JSON array"))
        end
    elseif c == '"'
        return _jstring(r)
    elseif c == 't'
        _jexpect(r, "true"); return true
    elseif c == 'f'
        _jexpect(r, "false"); return false
    elseif c == 'n'
        _jexpect(r, "null"); return nothing
    else
        start = r.i[]
        while r.i[] <= ncodeunits(r.s) && (r.s[r.i[]] in "+-.eE" || isdigit(r.s[r.i[]]))
            r.i[] += 1
        end
        tok = r.s[start:r.i[]-1]
        v = tryparse(Int, tok)
        v === nothing && (v = tryparse(Float64, tok))
        v === nothing && throw(ArgumentError("invalid JSON number $(tok)"))
        return v
    end
end

function _jstring(r)
    _jpeek(r) == '"' || throw(ArgumentError("expected string in JSON"))
    _jnext!(r)
    io = IOBuffer()
    while true
        c = _jnext!(r)
        c == '"' && break
        if c == '\\'
            e = _jnext!(r)
            if e == 'u'
                code = parse(UInt32, r.s[r.i[]:r.i[]+3]; base = 16)
                r.i[] += 4
                if 0xd800 <= code <= 0xdbff && r.s[r.i[]] == '\\'
                    low = parse(UInt32, r.s[r.i[]+2:r.i[]+5]; base = 16)
                    r.i[] += 6
                    code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00)
                end
                write(io, Char(code))
            else
                write(io, e == 'n' ? '\n' : e == 't' ? '\t' : e == 'r' ? '\r' : e == 'b' ? '\b' : e == 'f' ? '\f' : e)
            end
        else
            write(io, c)
        end
    end
    return String(take!(io))
end
