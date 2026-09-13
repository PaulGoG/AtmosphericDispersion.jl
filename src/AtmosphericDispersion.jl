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

using QuadGK: quadgk
using TOML: TOML

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
export WindProfileSurface,
    SURFACE_WATER,
    SURFACE_AGRICULTURAL,
    SURFACE_FOREST_URBAN,
    WIND_PROFILE_SURFACES,
    RoughnessClass,
    ROUGHNESS_GRASSLAND_WATER,
    ROUGHNESS_ARABLE,
    ROUGHNESS_PASTURE,
    ROUGHNESS_RURAL,
    ROUGHNESS_FOREST_URBAN,
    ROUGHNESS_METROPOLIS,
    ROUGHNESS_CLASSES,
    PrecipitationType,
    PRECIPITATION_RAIN,
    PRECIPITATION_SNOW,
    PRECIPITATION_TYPES
export VerticalShapeCoefficients,
    vertical_shape_coefficients,
    RoughnessCoefficients,
    roughness_coefficients,
    roughness_length,
    lateral_coefficient,
    profile_exponent,
    WashoutCoefficients,
    washout_coefficients,
    PRECIPITATION_RATES
export wind_speed,
    lateral_dispersion,
    vertical_shape,
    roughness_correction,
    vertical_dispersion,
    REFERENCE_HEIGHT,
    PROFILE_CEILING,
    SHORT_RELEASE_REFERENCE,
    SMOOTH_ROUGHNESS_LIMIT

export StackSource,
    Atmosphere,
    buoyancy_flux,
    momentum_flux,
    stability_parameter,
    STANDARD_GRAVITY,
    DRY_AIR_SPECIFIC_HEAT
export buoyancy_transition_distance,
    final_buoyant_rise,
    buoyant_rise,
    final_momentum_rise,
    momentum_rise,
    combined_rise,
    RiseCoefficients,
    BRIGGS_RISE,
    XOQDOQ_RISE,
    THESIS_RISE,
    plume_rise,
    BUOYANCY_FLUX_BREAKPOINT,
    MECHANISM_BALANCE_TOLERANCE

include("sectors.jl")
include("stability.jl")
include("windrose.jl")
include("surfaces.jl")
include("tables.jl")
include("dispersion.jl")
export Building,
    BuildingEnvelope,
    equivalent_height,
    equivalent_area,
    wake_broadened,
    WAKE_INFLUENCE_RADII,
    DEFAULT_WAKE_COEFFICIENT,
    NORMATIVE_WAKE_COEFFICIENT
export Site,
    mean_wind_speed,
    downwash_height,
    wake_height,
    release_height,
    transport_wind_speed,
    effective_height,
    corrected_lateral_dispersion,
    corrected_vertical_dispersion,
    WAKE_WIND_THRESHOLD
export plume_frame, dilution_instantaneous, dilution_extended, dilution_long_term
export DepositionVelocity, Nuclide, TRITIATED_WATER, TRITIUM_GAS, TRITIUM_DECAY_CONSTANT
export ResuspensionModel, RESUSPENSION_IAEA_SS57, RESUSPENSION_MAXWELL_ANSPAUGH
export WashoutModel, WASHOUT_NORMATIVE, WASHOUT_HTO, ogram_snow_washout
export WashoutSpecies, WASHOUT_TRITIUM_IODINE, WASHOUT_OTHER_NUCLIDES
export decay_factor,
    depletion_integral,
    dry_depletion_factor,
    wet_depletion_factor,
    depletion_factor,
    DEPLETION_INTEGRAL_FLOOR
export dry_deposition,
    wet_deposition,
    wet_deposition_sector,
    resuspension_factor,
    resuspended_concentration,
    RESUSPENSION_COEFFICIENTS
export ConfigurationError, RunConfiguration, load_configuration, configuration_from

include("source.jl")
include("plumerise.jl")
include("buildings.jl")
include("site.jl")
include("nuclides.jl")
include("depletion.jl")
include("dilution.jl")
include("deposition.jl")
include("config.jl")

end
