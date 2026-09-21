#=
Atmospheric dilution factors χ/Q, in s/m³.

Three regimes, distinguished by how much of the wind's variability the release
outlasts:

  instantaneous  minutes to an hour   the plume keeps a direction; a full
                                      three-dimensional Gaussian field
  extended       an hour to a day     the direction meanders within a sector;
                                      crosswind-averaged over that sector
  long term      beyond a day         the direction samples the whole rose;
                                      averaged over sectors and weighted by
                                      their frequencies

Only the long-term regime touches the wind rose, and so only it is exposed to
the direction convention.
=#

"""
    dilution_instantaneous(east, north, z, site, class, wind_bearing;
                           release_duration = SHORT_RELEASE_REFERENCE,
                           nuclide = nothing, washout = nothing)

Dilution factor χ/Q in s/m³ at the receptor `(east, north, z)` in metres, for a
release short enough that the wind holds a single direction `wind_bearing`,
given in radians clockwise from north as the direction the wind blows **from**.

The field is the Gaussian plume with total reflection at the ground,

    χ/Q = exp(−y²/2Σ_y²) [exp(−(z−H)²/2Σ_z²) + exp(−(z+H)²/2Σ_z²)] / (2π Σ_y Σ_z u)

with `x` the downwind and `y` the crosswind coordinate obtained by rotating the
receptor into the plume frame. Upwind of the source the factor is zero: this
model carries no upwind diffusion.
"""
function dilution_instantaneous(
    east::Real,
    north::Real,
    z::Real,
    site::AbstractSite,
    class::PasquillClass,
    wind_bearing::Real;
    release_duration::Real = SHORT_RELEASE_REFERENCE,
    nuclide::Union{Nothing,Nuclide} = nothing,
    washout::Union{Nothing,WashoutEvent} = nothing,
)
    z ≥ 0 || throw(DomainError(z, "receptor height cannot be below ground"))
    x, y = plume_frame(east, north, wind_bearing)
    x > 0 || return 0.0

    H = effective_height(x, site, class)
    Σy = corrected_lateral_dispersion(x, site, class; release_duration)
    Σz = corrected_vertical_dispersion(x, site, class)
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0

    D = _depletion(x, site, class, nuclide, washout)
    crosswind = exp(-y^2 / (2Σy^2)) / (sqrt(2π) * Σy)
    vertical = vertical_factor(z, H, Σz, mixing_layer(site), class)
    return D * crosswind * vertical / u
end

"""
    plume_frame(east, north, wind_bearing)

Downwind and crosswind coordinates `(x, y)` in metres of a receptor at
`(east, north)`, for a wind blowing **from** the bearing `wind_bearing`
(radians clockwise from north).

The plume axis points along the direction the wind blows towards, so a wind
from the north puts the axis due south. `x` is positive downwind of the source
and negative upwind of it.
"""
function plume_frame(east::Real, north::Real, wind_bearing::Real)
    # Unit vector along the transport direction, half a turn from the bearing
    # the wind blows from.
    ex, ey = -sin(wind_bearing), -cos(wind_bearing)
    x = east * ex + north * ey       # downwind
    y = -east * ey + north * ex      # crosswind, positive to the left of the axis
    return x, y
end

"""
    dilution_extended(east, north, site, class, wind_bearing, sectors = SectorGrid(16);
                      nuclide = nothing, washout = nothing)

Ground-level dilution factor χ/Q in s/m³ for a release long enough that the
wind direction meanders across a sector but short enough that it does not
sample the rose,

    χ/Q = √(2/π) exp(−H²/2Σ_z²) / (Σ_z u x θ_L)

uniform in the crosswind direction within the sector centred on the plume axis,
and zero outside it.

`wind_bearing` is the direction the wind blows **from**, radians clockwise from
north. `sectors` sets the sector width and defaults to the sixteen cardinal
sectors.
"""
function dilution_extended(
    east::Real,
    north::Real,
    site::AbstractSite,
    class::PasquillClass,
    wind_bearing::Real,
    sectors::SectorGrid = SectorGrid(16);
    nuclide::Union{Nothing,Nuclide} = nothing,
    washout::Union{Nothing,WashoutEvent} = nothing,
)
    x, y = plume_frame(east, north, wind_bearing)
    x > 0 || return 0.0
    θ_L = sector_width(sectors)
    abs(y / x) ≤ tan(θ_L / 2) || return 0.0

    H = effective_height(x, site, class)
    Σz = corrected_vertical_dispersion(x, site, class)
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0

    D = _depletion(x, site, class, nuclide, washout)
    # The crosswind-integrated vertical factor, spread over the arc the sector
    # subtends. Without a lid this is √(2/π)exp(−H²/2Σ_z²)/Σ_z exactly.
    vertical = crosswind_integrated_factor(H, Σz, mixing_layer(site), class)
    return D * vertical / (u * x * θ_L)
end

"""
    dilution_long_term(east, north, site, rose;
                       nuclide = nothing, washout = nothing)

Ground-level dilution factor χ/Q in s/m³ for a release long enough that the
wind direction samples the whole rose,

    χ_k/Q = √(2/π) F_k/(r θ_L) Σ_i [ F_ki D_i(r) exp(−H²/2Σ_z²) / (Σ_z ū) ]

summed over the Pasquill classes, with `F_k` the frequency of wind blowing
**towards** the receptor's sector `k` and `F_ki` the fraction of that time
spent in class `i`.

That is the form without a mixing lid. Under one, `√(2/π) exp(−H²/2Σ_z²)/Σ_z`
is replaced class by class with [`crosswind_integrated_factor`](@ref), exactly
as in [`dilution_extended`](@ref): the long-term factor is the extended one
summed over classes and weighted by the rose, and the two must agree on the
vertical profile.

`D_i` is the depletion factor of class `i`, which is one unless a `nuclide` is
given. It sits **inside** the class sum because it depends on the class through
both the transport speed and the vertical dispersion. The 2021 code instead
formed a separate depletion term as an unweighted sum of six exponentials, one
per class, and so returned up to six in the limit of no deposition at all,
where a surviving fraction must tend to one.

`rose` supplies both, in the correct sense whichever convention it was built
from — that is the point of [`WindRose`](@ref) storing blowing-towards
frequencies. The sector is the one containing the receptor's bearing.

The distance `r` is radial, as in the sector-averaged formulation this comes
from. The 2021 code instead projected the receptor onto the axis of its sector,
which shortens the distance by up to `1 − cos(θ_L/2)`, about 1.9 % for sixteen
sectors, and correspondingly inflates the dilution factor.
"""
function dilution_long_term(
    east::Real,
    north::Real,
    site::AbstractSite,
    rose::WindRose;
    nuclide::Union{Nothing,Nuclide} = nothing,
    washout::Union{Nothing,WashoutEvent} = nothing,
)
    r = hypot(east, north)
    r > 0 || return 0.0

    g = grid(rose)
    k = sector_of(g, east, north)
    F_k = frequency_toward(rose, k)
    iszero(F_k) && return 0.0

    total = 0.0
    for class in PASQUILL_CLASSES
        F_ki = stability_fraction(rose, k, class)
        iszero(F_ki) && continue
        H = effective_height(r, site, class)
        Σz = corrected_vertical_dispersion(r, site, class)
        u = transport_wind_speed(site, class)
        u > 0 || continue
        D = _depletion(r, site, class, nuclide, washout)
        # The same vertical factor as the extended regime, of which this is the
        # frequency-weighted sum, so the mixing lid enters here as it does there.
        vertical = crosswind_integrated_factor(H, Σz, mixing_layer(site), class)
        total += F_ki * D * vertical / u
    end

    return F_k * total / (r * sector_width(g))
end

# Dispatched rather than branched, so the undepleted path stays free of the
# nuclide machinery and both paths are type-stable. Washout without a nuclide is
# refused: the washout row is a property of the species.
function _depletion(::Real, ::AbstractSite, ::PasquillClass, ::Nothing, washout)
    washout === nothing ||
        throw(ArgumentError("a washout event needs a nuclide to act on; pass `nuclide`"))
    return 1.0
end

_depletion(
    x::Real,
    site::AbstractSite,
    class::PasquillClass,
    nuclide::Nuclide,
    washout::Union{Nothing,WashoutEvent},
) = depletion_factor(x, site, class, nuclide; washout)
