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
    AbstractSite

What the dilution, depletion and deposition functions ask of a release: its
[`effective_height`](@ref), [`transport_wind_speed`](@ref),
[`corrected_lateral_dispersion`](@ref), [`corrected_vertical_dispersion`](@ref)
and [`mixing_layer`](@ref). [`Site`](@ref) derives them from a stack and its
surroundings; [`PrescribedPlume`](@ref) takes them as given.
"""
abstract type AbstractSite end

"""
    Site(; source, atmosphere, buildings = BuildingEnvelope(), rise = BRIGGS_RISE,
         mixing = MIXING_TABULATED, dispersion = DISPERSION_HOSKER,
         stable_rise_wind = WIND_MEAN_OVER_RISE)

A stack in its surroundings: everything a dispersion calculation needs that
does not vary over the receptor grid.

`dispersion` selects the [`DispersionScheme`](@ref) the site's parameters are
evaluated with, and `stable_rise_wind` the [`StableRiseWind`](@ref) of its
stable final rise.

The release height, the buoyancy and momentum fluxes and the stability
parameter are all independent of receptor position and of stability class, so
they are computed once here. The 2021 code recomputed them inside the innermost
loop — the building reduction and the downwash correction, the latter six
wind-profile evaluations deep, were repeated for every point of every field.
"""
struct Site <: AbstractSite
    source::StackSource
    atmosphere::Atmosphere
    buildings::BuildingEnvelope
    rise::RiseCoefficients
    mixing::MixingLayer
    dispersion::DispersionScheme
    stable_rise_wind::StableRiseWind
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
            dispersion::DispersionScheme = DISPERSION_HOSKER,
            stable_rise_wind::StableRiseWind = WIND_MEAN_OVER_RISE,
    )
        return new(
            source,
            atmosphere,
            buildings,
            rise,
            mixing,
            dispersion,
            stable_rise_wind,
            wake_height(source, atmosphere, buildings),
            buoyancy_flux(source, atmosphere),
            momentum_flux(source, atmosphere),
            stability_parameter(atmosphere),
        )
    end
end

"""
    dispersion_scheme(site)

The [`DispersionScheme`](@ref) the site's parameters are evaluated with.
"""
dispersion_scheme(site::Site) = site.dispersion

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
    z = max(release_height(site), REFERENCE_HEIGHT)
    return wind_speed(site.atmosphere.reference_speed, z, site.atmosphere.surface, class)
end

"""
    STABLE_RISE_TOLERANCE

Relative change of the stable final rise between two iterations below which
[`stable_rise_wind_speed`](@ref) takes the rise and its mean wind as solved.
"""
const STABLE_RISE_TOLERANCE = 1e-12

"""
    STABLE_RISE_MAX_ITERATIONS

Iterations [`stable_rise_wind_speed`](@ref) allows the rise and its mean wind
to converge in before it gives up. The map contracts strongly — the rise
varies as the inverse cube root of the wind — and a handful suffice.
"""
const STABLE_RISE_MAX_ITERATIONS = 100

"""
    stable_rise_wind_speed(site, class)

The wind speed in m/s that the stable final rise of `site` is evaluated with in
the given stability class.

Under `WIND_MEAN_OVER_RISE` it is the mean of the wind profile between
the release height and the top of the stable rise, `Δh = c [F/(ūS)]^(1/3)`,
which depends on that mean in turn: the two are solved together by iterating
from the rise at the release-height wind, each pass averaging the profile over
the current rise and re-evaluating it, until the rise settles to
[`STABLE_RISE_TOLERANCE`](@ref). The profile has no skill below the height its
reference is quoted at, so the layer starts at [`REFERENCE_HEIGHT`](@ref) where
the release is lower, as [`transport_wind_speed`](@ref) is floored.

Under `WIND_AT_RELEASE_HEIGHT`, or where there is no stable rise to
average over — unstable or neutral air, or a plume without buoyancy — it is the
transport wind.
"""
function stable_rise_wind_speed(site::Site, class::PasquillClass)
    u = transport_wind_speed(site, class)
    site.stable_rise_wind == WIND_AT_RELEASE_HEIGHT && return u
    c, F, S = site.rise.stable_final, site.buoyancy, site.stability
    (S > 0 && F > 0 && c > 0 && u > 0) || return u
    atmosphere = site.atmosphere
    z₀ = max(release_height(site), REFERENCE_HEIGHT)
    Δh = c * (F / (u * S))^(1 / 3)
    for _ in 1:STABLE_RISE_MAX_ITERATIONS
        ū = layer_mean_wind_speed(
            atmosphere.reference_speed, z₀, z₀ + Δh, atmosphere.surface, class,)
        Δh′ = c * (F / (ū * S))^(1 / 3)
        abs(Δh′ - Δh) ≤ STABLE_RISE_TOLERANCE * Δh′ && return ū
        Δh = Δh′
    end
    error(
        "the stable final rise and the mean wind over it did not converge in " *
        "$STABLE_RISE_MAX_ITERATIONS iterations",
    )
end

"""
    plume_rise(x, site, class)

Rise of the plume above its release height, in metres, at downwind distance `x`
metres, using the transport wind of the given stability class, and for the
stable final rise the wind of [`stable_rise_wind_speed`](@ref).
"""
function plume_rise(x::Real, site::Site, class::PasquillClass)
    u = transport_wind_speed(site, class)
    stable_wind = stable_rise_wind_speed(site, class)
    w₀, D = site.source.exit_velocity, site.source.diameter
    return _plume_rise(
        x, site.buoyancy, site.momentum, site.stability, w₀, D, u, site.rise; stable_wind,)
end

"""
    effective_height(x, site, class)

Effective release height in metres at downwind distance `x` metres: the release
height after downwash and building wake, plus the plume rise attained by that
distance.
"""
effective_height(x::Real, site::Site, class::PasquillClass) = release_height(site) +
                                                              plume_rise(x, site, class)

"""
    corrected_lateral_dispersion(x, site, class; release_duration = SHORT_RELEASE_REFERENCE)

Lateral dispersion parameter in metres of the site's
[`DispersionScheme`](@ref), broadened by the building wake.
"""
function corrected_lateral_dispersion(
        x::Real,
        site::Site,
        class::PasquillClass;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    σy = lateral_dispersion(x, class, site.dispersion; release_duration)
    return wake_broadened(σy, effective_height(x, site, class), site.buildings)
end

"""
    corrected_vertical_dispersion(x, site, class)

Vertical dispersion parameter in metres of the site's
[`DispersionScheme`](@ref), broadened by the building wake.
"""
function corrected_vertical_dispersion(x::Real, site::Site, class::PasquillClass)
    σz = vertical_dispersion(x, class, site.dispersion, site.atmosphere.roughness)
    return wake_broadened(σz, effective_height(x, site, class), site.buildings)
end

"""
    mixing_layer(site)

The [`MixingLayer`](@ref) over the site.
"""
mixing_layer(site::Site) = site.mixing

"""
    PrescribedPlume(site; height = nothing, wind = nothing,
                    lateral = nothing, vertical = nothing)

The plume over `site` with some of its parameters stated rather than derived:
the effective height in metres, the transport wind speed in m/s, and the lateral
and vertical dispersion parameters in metres. Whatever is left as `nothing`
comes from `site` as usual, and so does the mixing layer.

Published benchmarks are posed this way. Turner's worked problems give `H`, `u`
and the two σ read off his figures; the HPA-RPD-058 depletion tables give `H`
and the wind speed at it. Those numbers are the inputs of the problem, and
reproducing it means taking them as such rather than finding a stack and a
roughness class that happen to return them.

A prescribed `vertical` is a statement about one distance, so
[`depletion_integral`](@ref), which needs `σ_z` along the whole path, rejects
it.
"""
struct PrescribedPlume <: AbstractSite
    site::Site
    height::Union{Nothing,Float64}
    wind::Union{Nothing,Float64}
    lateral::Union{Nothing,Float64}
    vertical::Union{Nothing,Float64}

    function PrescribedPlume(
            site::Site;
            height::Union{Nothing,Real} = nothing,
            wind::Union{Nothing,Real} = nothing,
            lateral::Union{Nothing,Real} = nothing,
            vertical::Union{Nothing,Real} = nothing,
    )
        height === nothing ||
            height ≥ 0 ||
            throw(ArgumentError("a prescribed effective height cannot be negative"))
        wind === nothing ||
            wind > 0 ||
            throw(ArgumentError("a prescribed transport wind speed must be positive"))
        lateral === nothing ||
            lateral > 0 ||
            throw(
                ArgumentError("a prescribed lateral dispersion parameter must be positive"),
            )
        vertical === nothing ||
            vertical > 0 ||
            throw(
                ArgumentError(
                "a prescribed vertical dispersion parameter must be positive",
            ),
            )
        _float(v) = v === nothing ? nothing : Float64(v)
        return new(site, _float(height), _float(wind), _float(lateral), _float(vertical))
    end
end

mixing_layer(plume::PrescribedPlume) = mixing_layer(plume.site)
dispersion_scheme(plume::PrescribedPlume) = dispersion_scheme(plume.site)

_has_prescribed_vertical(::Site) = false
_has_prescribed_vertical(plume::PrescribedPlume) = plume.vertical !== nothing

# Each field is bound to a local before the test so that `=== nothing` narrows
# the union and the return type stays Float64.
function transport_wind_speed(plume::PrescribedPlume, class::PasquillClass)
    wind = plume.wind
    return wind === nothing ? transport_wind_speed(plume.site, class) : wind
end

function effective_height(x::Real, plume::PrescribedPlume, class::PasquillClass)
    height = plume.height
    return height === nothing ? effective_height(x, plume.site, class) : height
end

function corrected_lateral_dispersion(
        x::Real,
        plume::PrescribedPlume,
        class::PasquillClass;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    lateral = plume.lateral
    lateral === nothing || return lateral
    # The wake test takes the plume's own height, which may itself be prescribed.
    σy = lateral_dispersion(x, class, plume.site.dispersion; release_duration)
    return wake_broadened(σy, effective_height(x, plume, class), plume.site.buildings)
end

function corrected_vertical_dispersion(
        x::Real,
        plume::PrescribedPlume,
        class::PasquillClass,
)
    vertical = plume.vertical
    vertical === nothing || return vertical
    site = plume.site
    σz = vertical_dispersion(x, class, site.dispersion, site.atmosphere.roughness)
    return wake_broadened(σz, effective_height(x, plume, class), site.buildings)
end
