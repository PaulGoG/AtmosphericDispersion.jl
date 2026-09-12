"""
    AtmosphericDispersion

Gaussian-plume dispersion of stack releases: sector geometry, wind-rose
statistics, plume rise, Briggs dispersion parameters, depletion and ground
deposition.

The public interface is organised around explicit input types rather than
loose numbers, because the quantities this model consumes carry conventions that
are invisible in their magnitudes — a wind direction read in the wrong sense, or
a sector index counted the wrong way round, produces a result that is wrong by a
rotation and right in every other respect.
"""
module AtmosphericDispersion

export SectorGrid,
    nsectors,
    sector_width,
    sector_of,
    sector_bearing,
    sector_bounds,
    sector_name,
    opposite,
    bearing,
    CARDINAL_16
export PasquillClass,
    PASQUILL_A,
    PASQUILL_B,
    PASQUILL_C,
    PASQUILL_D,
    PASQUILL_E,
    PASQUILL_F,
    PASQUILL_CLASSES,
    pasquill,
    letter,
    classindex
export WindDirectionConvention,
    BlowingFrom,
    BlowingToward,
    WindRose,
    frequency_toward,
    frequency_from,
    stability_fraction,
    grid

include("sectors.jl")
include("stability.jl")
include("windrose.jl")

end
