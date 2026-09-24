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
    layer_mean_wind_speed(u₁₀, z₁, z₂, surface, class)

Mean of the power-law wind profile over the layer `z₁ ≤ z ≤ z₂`, in m/s,

    ū = (1/(z₂ − z₁)) ∫ u(z) dz

integrated in closed form, with the profile held at its ceiling value above
[`PROFILE_CEILING`](@ref) as [`wind_speed`](@ref) holds it. With `z₂ = z₁` it
is the wind speed at that height. This is the `u` of the stable final rise
under `WIND_MEAN_OVER_RISE`.
"""
function layer_mean_wind_speed(
        u₁₀::Real,
        z₁::Real,
        z₂::Real,
        surface::WindProfileSurface,
        class::PasquillClass,
)
    z₁ > 0 || throw(
        DomainError(z₁, "the wind profile is defined above ground, z₁ must be positive"),
    )
    z₂ ≥ z₁ || throw(DomainError(z₂, "the top of the layer cannot lie below its base"))
    u₁₀ ≥ 0 || throw(DomainError(u₁₀, "the reference wind speed cannot be negative"))
    a, b = float(z₁), float(z₂)
    a == b && return wind_speed(u₁₀, a, surface, class)
    m = profile_exponent(surface, class)
    # ∫ (z/z_r)^m dz between two heights below the ceiling.
    below(lo, hi) = REFERENCE_HEIGHT / (m + 1) *
                    ((hi / REFERENCE_HEIGHT)^(m + 1) - (lo / REFERENCE_HEIGHT)^(m + 1))
    c = PROFILE_CEILING
    ceiling = (c / REFERENCE_HEIGHT)^m
    integral = if b ≤ c
        below(a, b)
    elseif a ≥ c
        ceiling * (b - a)
    else
        below(a, c) + ceiling * (b - c)
    end
    return u₁₀ * integral / (b - a)
end

"""
    meander_broadened(σy, release_duration)

`σ_y` in metres broadened for a release longer than
[`SHORT_RELEASE_REFERENCE`](@ref) by the factor `(t_R / 600)^0.2`, which
accounts for the additional meander of the wind direction over the release;
unchanged at or below ten minutes. `release_duration` is in seconds. The rule
is the normative's and applies whichever [`DispersionScheme`](@ref) supplied
`σ_y`.
"""
function meander_broadened(σy::Real, release_duration::Real)
    release_duration > 0 ||
        throw(DomainError(release_duration, "release duration must be positive"))
    release_duration ≤ SHORT_RELEASE_REFERENCE && return float(σy)
    return σy * (release_duration / SHORT_RELEASE_REFERENCE)^0.2
end

# Briggs' formulae in the exact arithmetic of their exponents: no `^` for the
# half-integer powers, so that the open-country σ_y reproduces the 2021 code to
# the bit.
function _briggs(c::BriggsCoefficients, x::Real)
    y = 1 + c.b * x
    c.p == -0.5 && return c.a * x / sqrt(y)
    c.p == 0.5 && return c.a * x * sqrt(y)
    c.p == -1.0 && return c.a * x / y
    c.p == 0.0 && return c.a * float(x)
    return c.a * x * y^c.p
end

"""
    briggs_lateral_dispersion(x, class, scheme)

Lateral dispersion parameter `σ_y` in metres at downwind distance `x` metres
from Briggs' interpolation formulae [Briggs1973](@cite), as printed in the
Handbook on Atmospheric Diffusion [Hanna1982](@cite) Table 4.5,

    σ_y = a x (1 + b x)^(−1/2)

with `b = 10⁻⁴` under `DISPERSION_BRIGGS_OPEN_COUNTRY` and `4 × 10⁻⁴` under
`DISPERSION_BRIGGS_URBAN`; see [`briggs_lateral_coefficients`](@ref). Briggs
quotes the formulae for `10² < x < 10⁴ m`, [`validity_range`](@ref).
"""
function briggs_lateral_dispersion(x::Real, class::PasquillClass, scheme::DispersionScheme)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    return _briggs(briggs_lateral_coefficients(class, scheme), x)
end

"""
    briggs_vertical_dispersion(x, class, scheme)

Vertical dispersion parameter `σ_z` in metres at downwind distance `x` metres
from Briggs' interpolation formulae, `σ_z = a x (1 + b x)^p`, under
`DISPERSION_BRIGGS_OPEN_COUNTRY` or `DISPERSION_BRIGGS_URBAN`; see
[`briggs_vertical_coefficients`](@ref). Unlike Hosker's, these carry no
roughness dependence: the Handbook notes that they "are independent of release
height and roughness".
"""
function briggs_vertical_dispersion(x::Real, class::PasquillClass, scheme::DispersionScheme)
    x > 0 || throw(
        DomainError(x, "the vertical dispersion parameter is defined downwind of the source"),
    )
    return _briggs(briggs_vertical_coefficients(class, scheme), x)
end

"""
    eimutis_konicek_lateral_dispersion(x, class)

Lateral dispersion parameter `σ_y` in metres at downwind distance `x` metres
from the fit of [EimutisKonicek1972](@citet) to the Pasquill–Gifford curves,
`σ_y = a x^0.9031`; see [`eimutis_konicek_lateral_coefficient`](@ref). This is
the `σ_y` of NRC XOQDOQ [Sagendorf1982](@cite) and PAVAN [Bander1982](@cite).
"""
function eimutis_konicek_lateral_dispersion(x::Real, class::PasquillClass)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    return eimutis_konicek_lateral_coefficient(class) * x^EIMUTIS_KONICEK_LATERAL_EXPONENT
end

"""
    eimutis_konicek_vertical_dispersion(x, class)

Vertical dispersion parameter `σ_z` in metres at downwind distance `x` metres
from the fit of [EimutisKonicek1972](@citet) to the Pasquill–Gifford curves,
`σ_z = a x^b + c` with one set of coefficients on each of three ranges of
distance; see [`eimutis_konicek_vertical_coefficients`](@ref).
"""
function eimutis_konicek_vertical_dispersion(x::Real, class::PasquillClass)
    x > 0 || throw(
        DomainError(x, "the vertical dispersion parameter is defined downwind of the source"),
    )
    c = eimutis_konicek_vertical_coefficients(class, x)
    return c.a * x^c.b + c.c
end

"""
    lateral_dispersion(x, class; release_duration = SHORT_RELEASE_REFERENCE)
    lateral_dispersion(x, class, scheme; release_duration = SHORT_RELEASE_REFERENCE)

Lateral dispersion parameter `σ_y` in metres at downwind distance `x` metres.
Without a scheme it is Briggs' open-country form, the `σ_y` of the normative,

    σ_y = c₃ x / √(1 + 10⁻⁴ x)

and with one it is that scheme's; see [`DispersionScheme`](@ref). Either is
broadened for a release longer than ten minutes by
[`meander_broadened`](@ref). `release_duration` is in seconds.
"""
function lateral_dispersion(
        x::Real,
        class::PasquillClass;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    return lateral_dispersion(x, class, DISPERSION_HOSKER; release_duration)
end

function lateral_dispersion(
        x::Real,
        class::PasquillClass,
        scheme::DispersionScheme;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    σy = if scheme == DISPERSION_EIMUTIS_KONICEK
        eimutis_konicek_lateral_dispersion(x, class)
    elseif scheme == DISPERSION_HOSKER
        briggs_lateral_dispersion(x, class, DISPERSION_BRIGGS_OPEN_COUNTRY)
    else
        briggs_lateral_dispersion(x, class, scheme)
    end
    return meander_broadened(σy, release_duration)
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

"""
    vertical_dispersion(x, class, scheme, roughness)

Vertical dispersion parameter `σ_z` in metres at downwind distance `x` metres
under the given [`DispersionScheme`](@ref). `roughness` enters only through
Hosker's `σ_z`, the parameters of the other three schemes carrying no
roughness dependence.
"""
function vertical_dispersion(
        x::Real,
        class::PasquillClass,
        scheme::DispersionScheme,
        roughness::RoughnessClass,
)
    scheme == DISPERSION_HOSKER && return vertical_dispersion(x, class, roughness)
    scheme == DISPERSION_EIMUTIS_KONICEK &&
        return eimutis_konicek_vertical_dispersion(x, class)
    return briggs_vertical_dispersion(x, class, scheme)
end

"""
    dispersion_parameters(x, class, scheme, roughness;
                          release_duration = SHORT_RELEASE_REFERENCE)

Both dispersion parameters at downwind distance `x` metres under the given
[`DispersionScheme`](@ref), as a named tuple `(σy, σz)` in metres:
[`lateral_dispersion`](@ref) and [`vertical_dispersion`](@ref) of that scheme,
`σ_y` broadened for a long release.
"""
function dispersion_parameters(
        x::Real,
        class::PasquillClass,
        scheme::DispersionScheme,
        roughness::RoughnessClass;
        release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    return (
        σy = lateral_dispersion(x, class, scheme; release_duration),
        σz = vertical_dispersion(x, class, scheme, roughness),
    )
end

"""
    validity_range(scheme)

The range of downwind distance in metres, `(lower, upper)`, over which the
parameters of a [`DispersionScheme`](@ref) were fitted, as the intersection of
what its two parameters state; the ends are included.

  - Briggs quotes his formulae, open country and urban alike, for
    `10² < x < 10⁴ m` (Handbook Table 4.5 [Hanna1982](@cite)):
    `DISPERSION_BRIGGS_OPEN_COUNTRY` and `DISPERSION_BRIGGS_URBAN` give
    `(100, 10_000)`.
  - `DISPERSION_HOSKER` pairs Briggs' `σ_y` with Hosker's `σ_z`, which was
    fitted to Smith's curves "presented graphically out to distances of 100 km"
    (HPA-RPD-058 [SmithSimmonds2009](@cite) §3.2.2.1). The pairing is bounded
    by its `σ_y`: `(100, 10_000)`.
  - The Eimutis–Konicek fit reproduces the Pasquill–Gifford curves, which
    [Turner1970](@citet) draws from 100 m to 100 km: `(100, 100_000)`.

The evaluation itself refuses no distance. The configuration loader applies
the policy of `[model] extrapolation` to the receptor grid, and
[`within_validity`](@ref) tests a single distance.
"""
function validity_range(scheme::DispersionScheme)
    scheme == DISPERSION_EIMUTIS_KONICEK && return (100.0, 100_000.0)
    return (100.0, 10_000.0)
end

"""
    within_validity(x, scheme)

Whether the downwind distance `x` metres lies inside the
[`validity_range`](@ref) of the scheme, ends included.
"""
function within_validity(x::Real, scheme::DispersionScheme)
    lower, upper = validity_range(scheme)
    return lower ≤ x ≤ upper
end
