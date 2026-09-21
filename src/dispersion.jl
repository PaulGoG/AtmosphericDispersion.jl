"""
    REFERENCE_HEIGHT

Height in metres at which the reference wind speed of the power-law profile is
quoted, by convention ten metres.
"""
const REFERENCE_HEIGHT = 10.0

"""
    PROFILE_CEILING

Height in metres above which the power-law wind profile is held constant.
"""
const PROFILE_CEILING = 200.0

"""
    SHORT_RELEASE_REFERENCE

Reference release duration in seconds, ten minutes, below which no meander
broadening is applied to the lateral dispersion parameter.
"""
const SHORT_RELEASE_REFERENCE = 600.0

"""
    SMOOTH_ROUGHNESS_LIMIT

Roughness length in metres separating the two tabulated forms of the roughness
correction factor.
"""
const SMOOTH_ROUGHNESS_LIMIT = 0.1

"""
    wind_speed(u₁₀, z, surface, class)

Wind speed at height `z` metres, from the reference speed `u₁₀` at ten metres,
through the power law `u(z) = u₁₀ (z/10)^m`.

The exponent `m` comes from [`profile_exponent`](@ref). Above
[`PROFILE_CEILING`](@ref) the profile is held at its ceiling value: the power
law is a surface-layer fit, and extrapolating it into the free atmosphere
overstates the shear.

`z` must be positive and `u₁₀` non-negative.
"""
function wind_speed(u₁₀::Real, z::Real, surface::WindProfileSurface, class::PasquillClass)
    z > 0 || throw(
        DomainError(z, "the wind profile is defined above ground, z must be positive"),
    )
    u₁₀ ≥ 0 || throw(DomainError(u₁₀, "the reference wind speed cannot be negative"))
    m = profile_exponent(surface, class)
    return u₁₀ * (min(float(z), PROFILE_CEILING) / REFERENCE_HEIGHT)^m
end

"""
    lateral_dispersion(x, class; release_duration = SHORT_RELEASE_REFERENCE)

Lateral dispersion parameter `σ_y` in metres at downwind distance `x` metres,

    σ_y = c₃ x / √(1 + 10⁻⁴ x)

broadened for releases longer than ten minutes by the factor
`(t_R / 600)^0.2`, which accounts for the additional meander of the wind
direction over the release.

`release_duration` is in seconds.
"""
function lateral_dispersion(
        x::Real,
        class::PasquillClass;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    release_duration > 0 ||
        throw(DomainError(release_duration, "release duration must be positive"))
    σy = lateral_coefficient(class) * x / sqrt(1 + 1e-4 * x)
    release_duration ≤ SHORT_RELEASE_REFERENCE && return σy
    return σy * (release_duration / SHORT_RELEASE_REFERENCE)^0.2
end

"""
    vertical_shape(x, class)

Shape function `g(x) = a₁ x^b₁ / (1 + a₂ x^b₂)` of the vertical dispersion
parameter, in metres, at downwind distance `x` metres.
"""
function vertical_shape(x::Real, class::PasquillClass)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    c = vertical_shape_coefficients(class)
    return c.a₁ * x^c.b₁ / (1 + c.a₂ * x^c.b₂)
end

"""
    roughness_correction(x, roughness)

Roughness correction factor `F(x)` multiplying the shape function of the
vertical dispersion parameter, at downwind distance `x` metres.

Two forms are tabulated, selected by the roughness length: smooth surfaces
(`z₀ ≤ 0.1 m`) take `F = ln(c₁ x^d₁ / (1 + c₂ x^d₂))`, rough ones take
`F = ln(c₁ x^d₁ (1 + 1/(c₂ x^d₂)))`.
"""
function roughness_correction(x::Real, roughness::RoughnessClass)
    x > 0 ||
        throw(DomainError(x, "the roughness correction is defined downwind of the source"))
    c = roughness_coefficients(roughness)
    if c.z₀ > SMOOTH_ROUGHNESS_LIMIT
        return log(c.c₁ * x^c.d₁ * (1 + inv(c.c₂ * x^c.d₂)))
    end
    return log(c.c₁ * x^c.d₁ / (1 + c.c₂ * x^c.d₂))
end

"""
    vertical_dispersion(x, class, roughness)

Vertical dispersion parameter `σ_z` in metres at downwind distance `x` metres,
as the product of the shape function and the roughness correction,
`σ_z = g(x) F(x)`.

Throws a `DomainError` if the parameterisation returns a non-positive value.
That happens only at distances far shorter than the model is meant for — below
about a tenth of a millimetre for the smoothest class — but a negative `σ_z`
would otherwise propagate silently into every quantity downstream.
"""
function vertical_dispersion(x::Real, class::PasquillClass, roughness::RoughnessClass)
    σz = vertical_shape(x, class) * roughness_correction(x, roughness)
    σz > 0 || throw(
        DomainError(
        x,
        "the vertical dispersion parameterisation gives σ_z = $σz m at this distance, " *
        "which is outside its range of validity",
    ),
    )
    return σz
end
