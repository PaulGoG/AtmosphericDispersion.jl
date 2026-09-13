#=
Depletion of the airborne plume.

Three processes remove material from the plume as it travels: radioactive
decay, dry deposition to the ground, and washout by precipitation. Each is
expressed as a factor in [0, 1] multiplying the undepleted concentration, and
they compose by multiplication — they act on the same material in sequence, so
the surviving fraction is the product of the surviving fractions.

That is the correction to the 2021 code. There, the wet and dry factors were
*added*:

    χ = χ/Q · Q · DEC · (DEP_w + DEP_d)

Adding survival fractions has no limit that makes sense. Switch the rain off
and switch the deposition velocity off, so that nothing is removed at all, and
each factor tends to one while their sum tends to two: the concentration comes
out twice the undepleted value. See `depletion_factor` for the multiplicative
form and the suite for the limit that pins it.
=#

"""
    DEPLETION_INTEGRAL_FLOOR

Lower limit in metres of the dry-depletion integral.

The integrand carries a `1/σ_z` that diverges as the source is approached, so
the integral is started a short distance downwind. The 2021 code used the same
one-metre floor.
"""
const DEPLETION_INTEGRAL_FLOOR = 1.0

"""
    decay_factor(x, u, nuclide)

Surviving fraction after radioactive decay over a travel distance `x` metres at
transport speed `u` m/s, `exp(−λ x / u)`.
"""
function decay_factor(x::Real, u::Real, nuclide::Nuclide)
    x ≥ 0 || throw(DomainError(x, "travel distance cannot be negative"))
    u > 0 || throw(DomainError(u, "transport speed must be positive"))
    return exp(-nuclide.decay_constant * x / u)
end

"""
    depletion_integral(x, site, class)

The integral

    ∫ √(2/π) exp(−H(x')² / 2Σ_z(x')²) / Σ_z(x') dx'

from [`DEPLETION_INTEGRAL_FLOOR`](@ref) to `x`, in units of m⁻¹ × m, which sets
how much of the plume has met the ground by distance `x`.

Under a mixing layer the integrand carries the lid's image terms as well, which
is HPA-RPD-058 Eq. (3.17); it is [`crosswind_integrated_factor`](@ref) evaluated
along the plume axis.

Evaluated by adaptive Gauss–Kronrod quadrature. The 2021 code used a fifty-point
trapezoid rule, which is not adaptive and resolves the near field — where the
integrand varies fastest — worst.
"""
function depletion_integral(x::Real, site::Site, class::PasquillClass)
    # Narrowed to a concrete type before the quadrature: `float` of an abstract
    # `Real` is uninferable, and the whole of QuadGK's generic machinery then
    # has to be resolved dynamically.
    upper = Float64(x)
    upper > DEPLETION_INTEGRAL_FLOOR || return 0.0
    A = mixing_depth(class, site.mixing)
    integrand = function (x′::Float64)
        Σz = corrected_vertical_dispersion(x′, site, class)
        H = effective_height(x′, site, class)
        # The same crosswind-integrated vertical factor the dilution kernel
        # uses, so the lid's image terms enter the depletion too. HPA-RPD-058
        # Eq. (3.17) carries them; without them the plume keeps depositing as
        # though it could go on diluting upwards for ever.
        return crosswind_integrated_factor(H, Σz, A)
    end
    value, _ = quadgk(integrand, DEPLETION_INTEGRAL_FLOOR, upper)
    return value
end

"""
    dry_depletion_factor(x, site, class, nuclide)

Surviving airborne fraction after dry deposition over `x` metres,

    exp[ −√(2/π) (v_d / u) ∫ exp(−H²/2Σ_z²)/Σ_z dx' ]

The *low* deposition velocity of the nuclide is used, which removes the least
material and so is conservative for the airborne concentration.
"""
function dry_depletion_factor(x::Real, site::Site, class::PasquillClass, nuclide::Nuclide)
    v_d = nuclide.deposition_velocity.low
    iszero(v_d) && return 1.0
    u = transport_wind_speed(site, class)
    u > 0 || throw(DomainError(u, "transport speed must be positive"))
    # The √(2/π) is inside the integrand now, in crosswind_integrated_factor.
    return exp(-(v_d / u) * depletion_integral(x, site, class))
end

"""
    wet_depletion_factor(duration, precipitation, rate)

Surviving airborne fraction after washout over a precipitation event of
`duration` seconds, `exp(−Λ t)`.

The *low* washout coefficient is used, for the same reason the low deposition
velocity is: it removes the least material and so is conservative for the
airborne concentration. The high coefficient drives the ground deposition
instead.
"""
function wet_depletion_factor(
    duration::Real,
    precipitation::PrecipitationType,
    rate::Real,
    model::WashoutModel = WASHOUT_NORMATIVE,
)
    duration ≥ 0 || throw(DomainError(duration, "washout duration cannot be negative"))
    return exp(-washout_coefficients(precipitation, rate, model).low * duration)
end

"""
    depletion_factor(x, site, class, nuclide; washout_duration = 0, precipitation, rate)

Total surviving airborne fraction over `x` metres: the product of the decay,
dry-deposition and washout factors.

They multiply because they act in sequence on the same material. With no
washout, no deposition and no decay each factor is one and so is the product,
which is the limit that fixes the composition — the 2021 code added the wet and
dry factors and so returned two in that limit.

Washout applies only where `washout_duration` is positive; with the default of
zero the plume is taken to travel in dry air.
"""
function depletion_factor(
    x::Real,
    site::Site,
    class::PasquillClass,
    nuclide::Nuclide;
    washout_duration::Real = 0.0,
    precipitation::PrecipitationType = PRECIPITATION_RAIN,
    rate::Real = first(PRECIPITATION_RATES),
    washout_model::WashoutModel = WASHOUT_NORMATIVE,
)
    u = transport_wind_speed(site, class)
    f = decay_factor(x, u, nuclide) * dry_depletion_factor(x, site, class, nuclide)
    washout_duration > 0 &&
        (f *= wet_depletion_factor(washout_duration, precipitation, rate, washout_model))
    return f
end
