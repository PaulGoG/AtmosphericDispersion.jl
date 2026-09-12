"""
    STANDARD_GRAVITY

Standard acceleration of gravity, 9.80665 m/s², the SI defined value.
"""
const STANDARD_GRAVITY = 9.80665

"""
    DRY_AIR_SPECIFIC_HEAT

Isobaric specific heat capacity of dry air, 1005 J/(kg·K).

This is the value the atmospheric stability parameter needs: the term `g/c_p`
in it is the dry adiabatic lapse rate, a property of the ambient air the plume
rises through, not of the gas leaving the stack.
"""
const DRY_AIR_SPECIFIC_HEAT = 1005.0

"""
    StackSource(; height, diameter, exit_velocity, exit_density, exit_temperature)

Geometry and discharge conditions of an elevated point source.

- `height` — stack height above ground, m
- `diameter` — inner diameter at the stack exit, m
- `exit_velocity` — vertical gas velocity at the exit, m/s
- `exit_density` — density of the emitted gas at the exit, kg/m³
- `exit_temperature` — temperature of the emitted gas, K

The buoyancy and momentum fluxes that drive plume rise are derived from these
together with the ambient state; see [`buoyancy_flux`](@ref) and
[`momentum_flux`](@ref).
"""
struct StackSource
    height::Float64
    diameter::Float64
    exit_velocity::Float64
    exit_density::Float64
    exit_temperature::Float64

    function StackSource(;
        height::Real,
        diameter::Real,
        exit_velocity::Real,
        exit_density::Real,
        exit_temperature::Real,
    )
        height > 0 || throw(ArgumentError("stack height must be positive, got $height m"))
        diameter > 0 ||
            throw(ArgumentError("stack diameter must be positive, got $diameter m"))
        exit_velocity ≥ 0 ||
            throw(ArgumentError("exit velocity cannot be negative, got $exit_velocity m/s"))
        exit_density > 0 ||
            throw(ArgumentError("exit density must be positive, got $exit_density kg/m³"))
        exit_temperature > 0 || throw(
            ArgumentError("exit temperature must be positive, got $exit_temperature K"),
        )
        return new(height, diameter, exit_velocity, exit_density, exit_temperature)
    end
end

"""
    Atmosphere(; reference_speed, temperature, density, lapse_rate, surface, roughness,
                 specific_heat = DRY_AIR_SPECIFIC_HEAT)

Ambient state the plume disperses into.

- `reference_speed` — wind speed at [`REFERENCE_HEIGHT`](@ref), m/s
- `temperature` — ambient air temperature, K
- `density` — ambient air density, kg/m³
- `lapse_rate` — vertical temperature gradient `dT/dz`, K/m, positive for an
  inversion
- `surface` — [`WindProfileSurface`](@ref) driving the wind profile exponent
- `roughness` — [`RoughnessClass`](@ref) driving the vertical dispersion
- `specific_heat` — isobaric specific heat of the ambient air, J/(kg·K)

!!! note "Lapse rate and stability class are not independent"
    The Pasquill class and the lapse rate describe the same stratification. A
    physically consistent calculation varies `lapse_rate` with the class —
    negative for the unstable classes A to C, near the adiabatic value for D,
    positive for E and F. Holding one fixed while sweeping the other, as the
    2021 code did with a single strongly stable gradient applied to every class,
    makes the unstable classes internally inconsistent.
"""
struct Atmosphere
    reference_speed::Float64
    temperature::Float64
    density::Float64
    lapse_rate::Float64
    specific_heat::Float64
    surface::WindProfileSurface
    roughness::RoughnessClass

    function Atmosphere(;
        reference_speed::Real,
        temperature::Real,
        density::Real,
        lapse_rate::Real,
        surface::WindProfileSurface,
        roughness::RoughnessClass,
        specific_heat::Real = DRY_AIR_SPECIFIC_HEAT,
    )
        reference_speed ≥ 0 || throw(
            ArgumentError(
                "reference wind speed cannot be negative, got $reference_speed m/s",
            ),
        )
        temperature > 0 ||
            throw(ArgumentError("ambient temperature must be positive, got $temperature K"))
        density > 0 ||
            throw(ArgumentError("ambient density must be positive, got $density kg/m³"))
        isfinite(lapse_rate) ||
            throw(ArgumentError("lapse rate must be finite, got $lapse_rate K/m"))
        specific_heat > 0 || throw(
            ArgumentError("specific heat must be positive, got $specific_heat J/(kg·K)"),
        )
        return new(
            reference_speed,
            temperature,
            density,
            lapse_rate,
            specific_heat,
            surface,
            roughness,
        )
    end
end

"""
    buoyancy_flux(source, atmosphere)

Buoyancy flux parameter `F` in m⁴/s³,

    F = (ρ − ρ₀)/ρ · g · w₀ · (D/2)²

from the density deficit of the emitted gas against the ambient air.

Negative for a plume denser than the ambient air, in which case the rise
correlations of this model do not apply — they describe buoyant plumes — and
[`buoyant_rise`](@ref) will reject it.
"""
function buoyancy_flux(source::StackSource, atmosphere::Atmosphere)
    Δρ = atmosphere.density - source.exit_density
    return Δρ / atmosphere.density *
           STANDARD_GRAVITY *
           source.exit_velocity *
           (source.diameter / 2)^2
end

"""
    momentum_flux(source, atmosphere)

Momentum flux parameter `F_m` in m⁴/s²,

    F_m = ρ₀/ρ · w₀² · (D/2)²
"""
function momentum_flux(source::StackSource, atmosphere::Atmosphere)
    return source.exit_density / atmosphere.density *
           source.exit_velocity^2 *
           (source.diameter / 2)^2
end

"""
    stability_parameter(atmosphere)

Atmospheric stability parameter `S` in s⁻²,

    S = g/T · (g/c_p + dT/dz)

the squared Brunt–Väisälä frequency of the ambient stratification. Positive in
stably stratified air, zero at the dry adiabatic lapse rate, negative in
unstable air.

The specific heat entering through `g/c_p` is that of the ambient air, which is
what makes that term the dry adiabatic lapse rate. The 2021 code used the
specific heat of the emitted water vapour here, 1871 J/(kg·K) instead of about
1005, which halves the adiabatic correction and shifts `S` by some fifteen per
cent for the gradient it ran with.
"""
function stability_parameter(atmosphere::Atmosphere)
    return STANDARD_GRAVITY / atmosphere.temperature *
           (STANDARD_GRAVITY / atmosphere.specific_heat + atmosphere.lapse_rate)
end
