"""
    WindDirectionConvention

Which of the two mutually inverse readings a tabulated wind direction carries.
Subtypes are [`BlowingFrom`](@ref) and [`BlowingToward`](@ref).

The distinction is a 180° rotation of the entire rose. Getting it wrong rotates
a long-term dispersion field by half a turn while leaving every magnitude, every
sum rule and every unit intact, so nothing downstream can detect the error. The
convention is therefore a required argument of the [`WindRose`](@ref)
constructor, never a default.
"""
abstract type WindDirectionConvention end

"""
    BlowingFrom()

The direction a wind is blowing *from*: the meteorological convention, and what
ADMS reads from the `PHI` field of a `.met` file. A westerly — 270° — blows from
the west towards the east. Measured wind roses are published this way.
"""
struct BlowingFrom <: WindDirectionConvention end

"""
    BlowingToward()

The direction a wind is blowing *towards*. This is the sense the long-term
sector-averaged dispersion formula needs, because its directional weight asks
how often the plume reaches the receptor's sector.
"""
struct BlowingToward <: WindDirectionConvention end

"""
    WindRose(grid, frequencies, convention; stability)

Joint frequency distribution of wind direction and atmospheric stability over
the sectors of `grid`.

`frequencies[k]` is the fraction of the time the wind lies in sector `k`, read
in the sense given by `convention` — either `BlowingFrom()` or
`BlowingToward()`. Frequencies must be non-negative and sum to one.

`stability` is the distribution over Pasquill classes *conditional on the
sector*: either a `6`-element vector, when the stability distribution is taken
to be the same in every sector, or an `nsectors(grid) × 6` matrix whose rows
each sum to one. Column `i` corresponds to `PASQUILL_CLASSES[i]`.

Whatever convention the input carries, a `WindRose` stores blowing-towards
frequencies internally, so the formulae that consume it cannot be fed the wrong
sense. Query either reading with [`frequency_toward`](@ref) and
[`frequency_from`](@ref).

# Examples

```jldoctest
julia> grid = SectorGrid(16);

julia> rose = WindRose(grid, fill(1/16, 16), BlowingFrom(); stability = fill(1/6, 6));

julia> frequency_toward(rose, 1) ≈ 1/16
true
```
"""
struct WindRose{T<:AbstractFloat}
    grid::SectorGrid
    toward::Vector{T}
    stability::Matrix{T}

    function WindRose(
            grid::SectorGrid,
            frequencies::AbstractVector{<:Real},
            convention::WindDirectionConvention;
            stability::AbstractVecOrMat{<:Real},
    )
        n = nsectors(grid)
        length(frequencies) == n || throw(
            DimensionMismatch(
            "the rose has $(length(frequencies)) directional frequencies but the grid has $n sectors",
        ),
        )

        T = float(promote_type(eltype(frequencies), eltype(stability)))
        f = collect(T, frequencies)
        _checkdistribution(f, "directional frequencies")

        s = _stabilitymatrix(T, stability, n)
        for k in 1:n
            _checkdistribution(view(s, k, :), "stability fractions of sector $k")
        end

        return new{T}(grid, _toward(convention, grid, f), s)
    end
end

# Storage is always blowing-towards; a blowing-from rose is rotated by half a turn.
_toward(::BlowingToward, ::SectorGrid, f::Vector) = f
_toward(::BlowingFrom, grid::SectorGrid, f::Vector) = [f[opposite(grid, k)]
                                                       for k in 1:nsectors(grid)]

function _stabilitymatrix(::Type{T}, stability::AbstractVector{<:Real}, n::Int) where {T}
    length(stability) == length(PASQUILL_CLASSES) || throw(
        DimensionMismatch(
        "a sector-independent stability distribution needs $(length(PASQUILL_CLASSES)) entries, got $(length(stability))",
    ),
    )
    return repeat(reshape(collect(T, stability), 1, :), n, 1)
end

function _stabilitymatrix(::Type{T}, stability::AbstractMatrix{<:Real}, n::Int) where {T}
    size(stability) == (n, length(PASQUILL_CLASSES)) || throw(
        DimensionMismatch(
        "the stability matrix is $(size(stability)) but should be ($n, $(length(PASQUILL_CLASSES)))",
    ),
    )
    return collect(T, stability)
end

const _DISTRIBUTION_ATOL = 1e-3

function _checkdistribution(p::AbstractArray{<:AbstractFloat}, what::AbstractString)
    all(isfinite, p) || throw(ArgumentError("the $what contain a non-finite value"))
    all(≥(0), p) || throw(ArgumentError("the $what contain a negative value"))
    total = sum(p)
    abs(total - 1) ≤ _DISTRIBUTION_ATOL || throw(
        ArgumentError(
        "the $what sum to $total, not to 1; supply a normalised distribution rather than raw counts",
    ),
    )
    return nothing
end

"""
    grid(rose)

The [`SectorGrid`](@ref) the rose is defined on.
"""
grid(rose::WindRose) = rose.grid

nsectors(rose::WindRose) = nsectors(rose.grid)

"""
    frequency_toward(rose, k)

Fraction of the time the wind blows *towards* sector `k` — the directional
weight the long-term sector-averaged dispersion formula applies to a receptor
lying in sector `k`.
"""
function frequency_toward(rose::WindRose, k::Integer)
    _checksector(rose.grid, k)
    return rose.toward[k]
end

"""
    frequency_from(rose, k)

Fraction of the time the wind blows *from* sector `k`, the reading a measured
wind rose and an ADMS `.met` file carry.
"""
frequency_from(rose::WindRose, k::Integer) = frequency_toward(rose, opposite(rose.grid, k))

"""
    stability_fraction(rose, k, class)

Fraction of the time Pasquill `class` prevails, conditional on the wind blowing
towards sector `k`.
"""
function stability_fraction(rose::WindRose, k::Integer, class::PasquillClass)
    _checksector(rose.grid, k)
    return rose.stability[k, classindex(class)]
end

function Base.show(io::IO, ::MIME"text/plain", rose::WindRose{T}) where {T}
    n = nsectors(rose)
    println(io, "WindRose{", T, "} over ", n, " sectors")
    k = argmax(rose.toward)
    print(
        io,
        "  most frequent transport towards ",
        sector_name(rose.grid, k),
        " (",
        round(100 * rose.toward[k]; digits = 2),
        "% of the time), ",
        "i.e. wind from ",
        sector_name(rose.grid, opposite(rose.grid, k)),
    )
    return nothing
end
