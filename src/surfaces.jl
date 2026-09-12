"""
    WindProfileSurface

Surface category controlling the exponent of the power-law wind profile
`u(z) = u₁₀ (z/10)^m`.

This is a coarser classification than [`RoughnessClass`](@ref), which controls
the vertical dispersion parameter. The two are tabulated separately and are not
interchangeable: a site is described by one of each.
"""
@enum WindProfileSurface::UInt8 begin
    SURFACE_WATER = 1
    SURFACE_AGRICULTURAL = 2
    SURFACE_FOREST_URBAN = 3
end

"""
    WIND_PROFILE_SURFACES

The three wind-profile surface categories, in table order.
"""
const WIND_PROFILE_SURFACES = (SURFACE_WATER, SURFACE_AGRICULTURAL, SURFACE_FOREST_URBAN)

"""
    RoughnessClass

Surface roughness category controlling the vertical dispersion parameter
`σ_z`, from open water through to dense urban fabric.

Each class carries a roughness length `z₀` (see [`roughness_length`](@ref)) and
a set of coefficients for the roughness correction factor of `σ_z`.
"""
@enum RoughnessClass::UInt8 begin
    ROUGHNESS_GRASSLAND_WATER = 1
    ROUGHNESS_ARABLE = 2
    ROUGHNESS_PASTURE = 3
    ROUGHNESS_RURAL = 4
    ROUGHNESS_FOREST_URBAN = 5
    ROUGHNESS_METROPOLIS = 6
end

"""
    ROUGHNESS_CLASSES

The six roughness classes, in order of increasing roughness length.
"""
const ROUGHNESS_CLASSES = (
    ROUGHNESS_GRASSLAND_WATER,
    ROUGHNESS_ARABLE,
    ROUGHNESS_PASTURE,
    ROUGHNESS_RURAL,
    ROUGHNESS_FOREST_URBAN,
    ROUGHNESS_METROPOLIS,
)

"""
    PrecipitationType

Form of precipitation scavenging the plume. Rain removes airborne material
between two and three orders of magnitude faster than snow at the same
intensity.
"""
@enum PrecipitationType::UInt8 begin
    PRECIPITATION_RAIN = 1
    PRECIPITATION_SNOW = 2
end

"""
    PRECIPITATION_TYPES

Both precipitation types, rain first.
"""
const PRECIPITATION_TYPES = (PRECIPITATION_RAIN, PRECIPITATION_SNOW)
