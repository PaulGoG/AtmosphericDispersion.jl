"""
    SectorGrid(n)

A wind rose divided into `n` equal angular sectors.

Sectors are indexed `1:n` **clockwise from north**, with sector 1 **centred** on
north. For the usual `n = 16` this reproduces the cardinal sectors: sector 1 is
N and spans bearings 348.75°–11.25°, sector 2 is NNE, and so on. This is the
convention of ADMS meteorological input and of the joint-frequency tables that
long-term dispersion calculations are built on.

`n` must be even, so that every sector has an opposite (see [`opposite`](@ref));
converting a wind rose between the "blowing from" and "blowing toward"
conventions is a rotation by `n ÷ 2` sectors.

Centring sector 1 on north rather than starting it there is not cosmetic. It
places the cardinal directions at sector centres, maximally far from a sector
boundary, so binning them cannot be decided by floating-point round-off.

# Examples

```jldoctest
julia> g = SectorGrid(16);

julia> sector_of(g, deg2rad(0.0)), sector_of(g, deg2rad(90.0))
(1, 5)

julia> rad2deg(sector_bearing(g, 5))
90.0
```
"""
struct SectorGrid
    n::Int

    function SectorGrid(n::Integer)
        n ≥ 4 || throw(ArgumentError("a wind rose needs at least 4 sectors, got $n"))
        iseven(n) || throw(
            ArgumentError(
            "the number of sectors must be even so that every sector has an opposite, got $n",
        ),
        )
        return new(Int(n))
    end
end

"""
    nsectors(grid)

Number of sectors in `grid`.
"""
nsectors(grid::SectorGrid) = grid.n

"""
    sector_width(grid)

Angular width of one sector, in radians.
"""
sector_width(grid::SectorGrid) = 2π / grid.n

"""
    bearing(east, north)

Compass bearing of the point `(east, north)` as seen from the origin, in radians
clockwise from north, wrapped to `[0, 2π)`.

The arguments are ground-plane displacements in the geographic frame, not the
plume-aligned frame: `east` is the displacement towards the east, `north`
towards the north. The origin returns a bearing of zero.
"""
function bearing(east::Real, north::Real)
    iszero(east) &&
        iszero(north) &&
        return zero(float(promote_type(typeof(east), typeof(north))))
    return mod2pi(atan(east, north))
end

"""
    sector_of(grid, β)
    sector_of(grid, east, north)

Index of the sector of `grid` containing the bearing `β` (radians clockwise from
north), or containing the point `(east, north)`.

Binning rounds to the nearest sector centre rather than truncating a ratio.
A bearing at a sector centre therefore bins exactly, because the ratio is an
integer and a round-off of an ulp either way still rounds to it. This is what
makes the cardinal directions safe: they lie at the centres.

On a sector boundary the ratio is a half-integer, the tie rule takes the sector
of higher bearing, and an ulp of round-off can place it either side. That is
unavoidable — a boundary bearing is generally not representable — and harmless,
since the two sectors meeting there are equally defensible. Do not build a
result on which side a boundary falls.
"""
function sector_of(grid::SectorGrid, β::Real)
    k = round(Int, β / sector_width(grid), RoundNearestTiesUp)
    return mod(k, grid.n) + 1
end

sector_of(grid::SectorGrid, east::Real, north::Real) = sector_of(grid, bearing(east, north))

"""
    sector_bearing(grid, k)

Bearing of the centre of sector `k`, in radians clockwise from north.
"""
function sector_bearing(grid::SectorGrid, k::Integer)
    _checksector(grid, k)
    return (k - 1) * sector_width(grid)
end

"""
    sector_bounds(grid, k)

Bearings of the two edges of sector `k` as a tuple `(lower, upper)`, in radians
clockwise from north. For sector 1 the lower edge wraps, so `lower > upper`.
"""
function sector_bounds(grid::SectorGrid, k::Integer)
    c = sector_bearing(grid, k)
    h = sector_width(grid) / 2
    return (mod2pi(c - h), mod2pi(c + h))
end

"""
    opposite(grid, k)

Index of the sector diametrically opposite `k`.

Wind blowing *from* sector `k` blows *towards* `opposite(grid, k)`, which is
what makes this the conversion between the two wind-direction conventions.
"""
function opposite(grid::SectorGrid, k::Integer)
    _checksector(grid, k)
    return mod(k - 1 + grid.n ÷ 2, grid.n) + 1
end

function _checksector(grid::SectorGrid, k::Integer)
    1 ≤ k ≤ grid.n || throw(BoundsError("sector index $k outside 1:$(grid.n)"))
    return nothing
end

"""
    CARDINAL_16

Names of the sixteen cardinal sectors, in the order `SectorGrid(16)` indexes
them: clockwise from north.
"""
const CARDINAL_16 = (
    "N",
    "NNE",
    "NE",
    "ENE",
    "E",
    "ESE",
    "SE",
    "SSE",
    "S",
    "SSW",
    "SW",
    "WSW",
    "W",
    "WNW",
    "NW",
    "NNW",
)

"""
    sector_name(grid, k)

Cardinal name of sector `k` when `grid` has sixteen sectors, otherwise a bearing
label such as `"123.8°"`.
"""
function sector_name(grid::SectorGrid, k::Integer)
    _checksector(grid, k)
    grid.n == 16 && return CARDINAL_16[k]
    return string(round(rad2deg(sector_bearing(grid, k)); digits = 1), "°")
end
