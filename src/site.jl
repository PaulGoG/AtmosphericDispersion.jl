"""
    WAKE_WIND_THRESHOLD

Wind speed in m/s below which a plume clearing the buildings is taken to escape
their wake rather than be drawn down into it.
"""
const WAKE_WIND_THRESHOLD = 5.0

"""
    mean_wind_speed(atmosphere, z)

Wind speed at height `z` metres averaged over the six Pasquill classes.

The stack-downwash and building-wake corrections are geometric: they ask
whether the efflux can overcome the wind at the stack top at all, not what the
stratification is doing. They therefore take this class-averaged speed rather
than a class-specific one.
"""
function mean_wind_speed(atmosphere::Atmosphere, z::Real)
    total = 0.0
    for class in PASQUILL_CLASSES
        total += wind_speed(atmosphere.reference_speed, z, atmosphere.surface, class)
    end
    return total / length(PASQUILL_CLASSES)
end

"""
    downwash_height(source, atmosphere)

Stack height in metres corrected for aerodynamic downwash at the stack itself.

Where the efflux velocity fails to clear about one and a half times the wind
speed, the plume is drawn into the lee of the stack and the effective release
height drops by `2 (1.5 − w₀/u) D`.
"""
function downwash_height(source::StackSource, atmosphere::Atmosphere)
    u = mean_wind_speed(atmosphere, source.height)
    u > 0 || return source.height
    w₀ = source.exit_velocity
    w₀ < 1.5u || return source.height
    return source.height - 2 * (1.5 - w₀ / u) * source.diameter
end

"""
    wake_height(source, atmosphere, envelope)

Release height in metres after both the stack downwash and entrainment into the
aerodynamic cavity of nearby buildings.

A plume leaving below the building tops is taken into the cavity and released
at ground level. One clearing two and a half times the building height escapes
untouched, as does one in wind below [`WAKE_WIND_THRESHOLD`](@ref), too weak to
drive the entrainment. Between those the release height is reduced by
`1.5 H_b − 0.6 H₁`.
"""
function wake_height(
    source::StackSource,
    atmosphere::Atmosphere,
    envelope::BuildingEnvelope,
)
    H₁ = downwash_height(source, atmosphere)
    h = equivalent_height(envelope)
    H₁ < h && return 0.0
    H₁ > 2.5h && return H₁
    mean_wind_speed(atmosphere, H₁) < WAKE_WIND_THRESHOLD && return H₁
    return H₁ - (1.5h - 0.6H₁)
end

"""
    Site(; source, atmosphere, buildings = BuildingEnvelope(), rise, mixing,
         fixed_height = nothing, fixed_wind = nothing,
         fixed_lateral = nothing, fixed_vertical = nothing)

A stack in its surroundings: everything a dispersion calculation needs that
does not vary over the receptor grid.

The release height, the buoyancy and momentum fluxes and the stability
parameter are all independent of receptor position and of stability class, so
they are computed once here. The 2021 code recomputed them inside the innermost
loop — the building reduction and the downwash correction, the latter six
wind-profile evaluations deep, were repeated for every point of every field.

`fixed_height`, `fixed_wind`, `fixed_lateral` and `fixed_vertical` override the
effective release height in metres, the transport wind speed in m/s, and the two
dispersion parameters in metres — bypassing plume rise, the wind profile, the
dispersion parameterisation and the building wake.
Published benchmarks state both as inputs rather than deriving them — Turner's
worked problems, the HPA-RPD-058 depletion tables, the NRC test cases — so
reproducing one requires setting them, not reverse-engineering a stack geometry
and a roughness class that happen to produce them. Turner's worked problems
print σ_y and σ_z read off his figures, and a σ read off a graph is an input to
the problem, not something to be recomputed.
"""

struct Site
    source::StackSource
    atmosphere::Atmosphere
    buildings::BuildingEnvelope
    rise::RiseCoefficients
    mixing::MixingLayer
    fixed_height::Union{Nothing,Float64}
    fixed_wind::Union{Nothing,Float64}
    fixed_lateral::Union{Nothing,Float64}
    fixed_vertical::Union{Nothing,Float64}
    release_height::Float64
    buoyancy::Float64
    momentum::Float64
    stability::Float64

    function Site(;
        source::StackSource,
        atmosphere::Atmosphere,
        buildings::BuildingEnvelope = BuildingEnvelope(),
        rise::RiseCoefficients = BRIGGS_RISE,
        mixing::MixingLayer = MIXING_TABULATED,
        fixed_height::Union{Nothing,Real} = nothing,
        fixed_wind::Union{Nothing,Real} = nothing,
        fixed_lateral::Union{Nothing,Real} = nothing,
        fixed_vertical::Union{Nothing,Real} = nothing,
    )
        fixed_height === nothing ||
            fixed_height ≥ 0 ||
            throw(ArgumentError("a fixed release height cannot be negative"))
        fixed_wind === nothing ||
            fixed_wind > 0 ||
            throw(ArgumentError("a fixed transport wind speed must be positive"))
        fixed_lateral === nothing ||
            fixed_lateral > 0 ||
            throw(ArgumentError("a fixed lateral dispersion parameter must be positive"))
        fixed_vertical === nothing ||
            fixed_vertical > 0 ||
            throw(ArgumentError("a fixed vertical dispersion parameter must be positive"))
        return new(
            source,
            atmosphere,
            buildings,
            rise,
            mixing,
            fixed_height === nothing ? nothing : Float64(fixed_height),
            fixed_wind === nothing ? nothing : Float64(fixed_wind),
            fixed_lateral === nothing ? nothing : Float64(fixed_lateral),
            fixed_vertical === nothing ? nothing : Float64(fixed_vertical),
            wake_height(source, atmosphere, buildings),
            buoyancy_flux(source, atmosphere),
            momentum_flux(source, atmosphere),
            stability_parameter(atmosphere),
        )
    end
end

"""
    release_height(site)

Height in metres at which the plume is released, after downwash and building
wake but before plume rise. Zero where the plume is trapped in a building
cavity.
"""
release_height(site::Site) = site.release_height

"""
    transport_wind_speed(site, class)

Wind speed in m/s carrying the plume, evaluated at the release height.

The power-law profile has no skill below the height its reference is quoted at,
and returns zero at the ground, so the evaluation height is floored at
[`REFERENCE_HEIGHT`](@ref). Without that floor a plume trapped in a building
cavity — release height zero — would be transported at zero wind speed, and
every dilution factor divides by it.
"""
function transport_wind_speed(site::Site, class::PasquillClass)
    # Bound to a local so the `=== nothing` test narrows the union; returning
    # the field directly leaves the return type Union{Nothing,Float64}.
    fixed = site.fixed_wind
    fixed === nothing || return fixed
    z = max(release_height(site), REFERENCE_HEIGHT)
    return wind_speed(site.atmosphere.reference_speed, z, site.atmosphere.surface, class)
end

"""
    plume_rise(x, site, class)

Rise of the plume above its release height, in metres, at downwind distance `x`
metres, using the transport wind of the given stability class.
"""
function plume_rise(x::Real, site::Site, class::PasquillClass)
    u = transport_wind_speed(site, class)
    F, Fₘ, S = site.buoyancy, site.momentum, site.stability
    w₀, D = site.source.exit_velocity, site.source.diameter

    r = site.rise
    buoyant = final_buoyant_rise(F, u, S, r)
    momentum = final_momentum_rise(Fₘ, w₀, D, u, S, r)
    total = buoyant + momentum
    balanced =
        iszero(total) || 2 * abs(buoyant - momentum) / total ≤ MECHANISM_BALANCE_TOLERANCE
    balanced && return combined_rise(x, F, Fₘ, w₀, u, S, D, r)
    momentum > buoyant && return momentum_rise(x, Fₘ, w₀, D, u, S, r)
    return buoyant_rise(x, F, u, S, r)
end

"""
    effective_height(x, site, class)

Effective release height in metres at downwind distance `x` metres: the release
height after downwash and building wake, plus the plume rise attained by that
distance.
"""
function effective_height(x::Real, site::Site, class::PasquillClass)
    fixed = site.fixed_height
    fixed === nothing && return release_height(site) + plume_rise(x, site, class)
    return fixed
end

"""
    corrected_lateral_dispersion(x, site, class; release_duration = SHORT_RELEASE_REFERENCE)

Lateral dispersion parameter in metres, broadened by the building wake.
"""
function corrected_lateral_dispersion(
    x::Real,
    site::Site,
    class::PasquillClass;
    release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    fixed = site.fixed_lateral
    fixed === nothing || return fixed
    σy = lateral_dispersion(x, class; release_duration)
    return wake_broadened(σy, effective_height(x, site, class), site.buildings)
end

"""
    corrected_vertical_dispersion(x, site, class)

Vertical dispersion parameter in metres, broadened by the building wake.
"""
function corrected_vertical_dispersion(x::Real, site::Site, class::PasquillClass)
    fixed = site.fixed_vertical
    fixed === nothing || return fixed
    σz = vertical_dispersion(x, class, site.atmosphere.roughness)
    return wake_broadened(σz, effective_height(x, site, class), site.buildings)
end
