"""
    PasquillClass

Pasquill–Gifford atmospheric stability class, `PASQUILL_A` (very unstable)
through `PASQUILL_F` (moderately stable).

The dispersion parameters, the wind-speed profile exponent and the plume-rise
correlations are all keyed by this class.
"""
@enum PasquillClass::UInt8 begin
    PASQUILL_A = 1
    PASQUILL_B = 2
    PASQUILL_C = 3
    PASQUILL_D = 4
    PASQUILL_E = 5
    PASQUILL_F = 6
end

"""
    PASQUILL_CLASSES

The six Pasquill classes in order, `PASQUILL_A` through `PASQUILL_F`. Sums over
stability classes iterate this tuple.
"""
const PASQUILL_CLASSES = (
    PASQUILL_A, PASQUILL_B, PASQUILL_C, PASQUILL_D, PASQUILL_E, PASQUILL_F,)

"""
    pasquill(c)

The Pasquill class named by the letter `c`, given as a `Char` or a
single-character `AbstractString`, case-insensitively.

Tabulated coefficients are conventionally published against the letters, so
parsing them is part of reading any input table.

# Examples

```jldoctest
julia> pasquill('D') === PASQUILL_D
true

julia> pasquill("f") === PASQUILL_F
true
```
"""
function pasquill(c::Char)
    i = Int(uppercase(c)) - Int('A') + 1
    1 ≤ i ≤ 6 || throw(ArgumentError("no Pasquill class '$c'; expected one of A-F"))
    return PASQUILL_CLASSES[i]
end

function pasquill(s::AbstractString)
    length(s) == 1 ||
        throw(ArgumentError("a Pasquill class is a single letter A-F, got $(repr(s))"))
    return pasquill(first(s))
end

pasquill(c::PasquillClass) = c

"""
    letter(class)

The letter naming a [`PasquillClass`](@ref).
"""
letter(class::PasquillClass) = Char(Int('A') + Int(class) - 1)

"""
    classindex(class)

Position of `class` in [`PASQUILL_CLASSES`](@ref), `1` for `PASQUILL_A` through
`6` for `PASQUILL_F`. This is the index into any table keyed by stability class.
"""
classindex(class::PasquillClass) = Int(class)
