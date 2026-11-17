#=
Depletion of the airborne plume.

Three processes remove material from the plume as it travels: radioactive
decay, dry deposition to the ground, and washout by precipitation. Each is
expressed as a factor in [0, 1] multiplying the undepleted concentration, and
they compose by multiplication: they act on the same material in sequence, so
the surviving fraction is the product of the surviving fractions.

Multiplication corrects the 2021 code, which added the wet and dry factors:

    χ = χ/Q · Q · DEC · (DEP_w + DEP_d)

Adding survival fractions fails in the limit of no removal: with no rain and a
zero deposition velocity each factor tends to one while their sum tends to two,
and the concentration is twice the undepleted value. See `depletion_factor` for
the multiplicative form and the test suite for the check of that limit.
=#

"""
    DEPLETION_INTEGRAL_FLOOR

Lower limit in metres of the dry-depletion integral.

The integrand carries a `1/σ_z` that diverges as the source is approached, so
the integral is started a short distance downwind. The 2021 code used the same
one-metre floor. A numerical control rather than part of the model, so it is
not exported.
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
is HPA-RPD-058 [SmithSimmonds2009](@cite) Eq. (3.17); it is [`crosswind_integrated_factor`](@ref) evaluated
along the plume axis.

Evaluated by adaptive Gauss–Kronrod quadrature. The 2021 code used a fifty-point
trapezoid rule, which is not adaptive and is least accurate in the near field,
where the integrand varies fastest.
"""
function depletion_integral(x::Real, site::AbstractSite, class::PasquillClass)
    _has_prescribed_vertical(site) && throw(
        ArgumentError(
        "the depletion integral needs σ_z along the whole path; a prescribed σ_z holds at one distance",
    ),
    )
    # Narrowed to a concrete type before the quadrature: `float` of an abstract
    # `Real` is uninferable, and the whole of QuadGK's generic machinery then
    # has to be resolved dynamically.
    upper = Float64(x)
    upper > DEPLETION_INTEGRAL_FLOOR || return 0.0
    layer = mixing_layer(site)
    integrand = function (x′::Float64)
        Σz = corrected_vertical_dispersion(x′, site, class)
        H = effective_height(x′, site, class)
        # The same crosswind-integrated vertical factor the dilution kernel
        # uses, so the lid's image terms enter the depletion too. HPA-RPD-058
        # Eq. (3.17) includes them; without them the deposition is computed as
        # though the plume could dilute upwards without limit.
        return crosswind_integrated_factor(H, Σz, layer, class)
    end
    value, _ = quadgk(integrand, DEPLETION_INTEGRAL_FLOOR, upper)
    return value
end

"""
    dry_depletion_factor(x, site, class, nuclide)

Surviving airborne fraction after dry deposition over `x` metres,

    exp[ −√(2/π) (v_d / u) ∫ exp(−H²/2Σ_z²)/Σ_z dx' ]

The low deposition velocity of the nuclide is used, which removes the least
material and so is conservative for the airborne concentration.
"""
function dry_depletion_factor(
        x::Real,
        site::AbstractSite,
        class::PasquillClass,
        nuclide::Nuclide,
)
    v_d = nuclide.deposition_velocity.low
    iszero(v_d) && return 1.0
    u = transport_wind_speed(site, class)
    u > 0 || throw(DomainError(u, "transport speed must be positive"))
    # The √(2/π) is inside the integrand, in crosswind_integrated_factor.
    return exp(-(v_d / u) * depletion_integral(x, site, class))
end

"""
    WashoutEvent(; duration, release_duration = duration,
                 precipitation = PRECIPITATION_RAIN,
                 rate = first(PRECIPITATION_RATES), model = WASHOUT_NORMATIVE)

Precipitation scavenging the plume: its `precipitation` type, its intensity
`rate` in mm/h, its `duration` in seconds during a release of
`release_duration` seconds, and the [`WashoutModel`](@ref) its coefficients
are taken from. The intensity must be one of [`PRECIPITATION_RATES`](@ref), the
washout table being defined only there.

The precipitation is taken to cover the whole path of the material released
while it falls, the fraction `min(duration/release_duration, 1)` of the
release: that material is scavenged for its whole travel time to the receptor,
the rest not at all. For the airborne plume this is the limit of a
precipitation long compared with the travel time, whose range
[`wet_depletion_factor`](@ref) states; [`wet_deposition`](@ref) takes the bound
for any timing instead. The default `release_duration` makes the event last the
whole release. Both durations must be finite.

Which row of the table applies is a property of the released species, so it is
stored in the [`Nuclide`](@ref).
"""
struct WashoutEvent
    precipitation::PrecipitationType
    rate::Float64
    duration::Float64
    fraction::Float64
    model::WashoutModel

    function WashoutEvent(;
            duration::Real,
            release_duration::Real = duration,
            precipitation::PrecipitationType = PRECIPITATION_RAIN,
            rate::Real = first(PRECIPITATION_RATES),
            model::WashoutModel = WASHOUT_NORMATIVE,
    )
        0 ≤ duration < Inf || throw(
            ArgumentError("washout duration must be finite and non-negative, got $duration s"),
        )
        0 ≤ release_duration < Inf || throw(
            ArgumentError(
            "release duration must be finite and non-negative, got $release_duration s",
        ),
        )
        duration > 0 && release_duration == 0 &&
            throw(
                ArgumentError("a washout event of positive duration needs a release duration"),
            )
        float(rate) in PRECIPITATION_RATES || throw(
            ArgumentError(
            "the washout table is defined at intensities $(PRECIPITATION_RATES) mm/h, got $rate",
        ),
        )
        fraction = duration == 0 ? 0.0 : min(Float64(duration) / release_duration, 1.0)
        return new(precipitation, Float64(rate), Float64(duration), fraction, model)
    end
end

"""
    precipitation_fraction(event)

Fraction of the release that falls in the precipitation of the
[`WashoutEvent`](@ref), `min(duration/release_duration, 1)`.
"""
precipitation_fraction(event::WashoutEvent) = event.fraction

"""
    washout_coefficients(event, species = WASHOUT_TRITIUM_IODINE)

The [`WashoutCoefficients`](@ref) of a [`WashoutEvent`](@ref) for the given
[`WashoutSpecies`](@ref).
"""
washout_coefficients(
    event::WashoutEvent,
    species::WashoutSpecies = WASHOUT_TRITIUM_IODINE,
) = washout_coefficients(event.precipitation, event.rate, event.model, species)

"""
    wet_depletion_factor(x, u, event, species = WASHOUT_TRITIUM_IODINE)
    wet_depletion_factor(event, species = WASHOUT_TRITIUM_IODINE)

Surviving airborne fraction of the release at downwind distance `x` metres,
transported at `u` m/s, under the washout `event`,

    1 − f + f exp(−Λ x/u)

with `f` the [`precipitation_fraction`](@ref) of the event. The material
released during the precipitation has been scavenged for its travel time
`x/u`, which is the time of CNCAN NSR-23 [CNCAN2004](@cite) Eq. (35), "the
period of time during which the radioactive material is subjected to the
process", and of HPA-RPD-058 [SmithSimmonds2009](@cite) Eq. (3.13); the rest of
the release travels in dry air.

The form is the limit of a precipitation that covers the path for much longer
than the travel time. Of the precipitation fixed over the site that holds the
share `f`, the one that begins with the release leaves the most airborne, or,
if it outlasts the release, the one that ends with it: each parcel is
scavenged only until it arrives or the precipitation stops. For
`x/u ≤ t_w ≤ T_R` the scavenged share then survives better than `exp(−Λ x/u)`
by the relative amount `((exp(Λ x/u) − 1)/Λ − x/u)/t_w`, which is small only
while `t_w ≫ (exp(Λ x/u) − 1)/Λ`, or `t_w ≫ x/u` where `Λ x/u ≪ 1`. A
precipitation within the release that begins later, while earlier material is
airborne, removes more. For `Λ = 10⁻³ s⁻¹`, 900 s of
precipitation in a 3600 s release and `x/u = 3333 s`, the surviving fraction is
0.92 if the precipitation falls at the start of the release and 0.47 if at the
end, against 0.76 from this form. Omitting washout bounds the airborne
concentration for any timing, as EPA ISC3 [EPA1995](@cite) §1.4 and IAEA
SRS-19 [IAEA2001](@cite) §3.8 note. [`load_configuration`](@ref) warns when,
at the edge of the receptor grid, the precipitation timed to leave the most
airborne would leave more than 10 % more than this form.

The low washout coefficient is used, for the same reason the low deposition
velocity is: it removes the least material and so is conservative for the
airborne concentration. The high coefficient drives the ground deposition
instead; see [`wet_deposition`](@ref).

The two-argument form is `exp(−Λ t_w)`, the whole event applied to the whole
release at every distance, as the 2021 code applied it; it is kept to
reproduce those results.
"""
function wet_depletion_factor(
        x::Real,
        u::Real,
        event::WashoutEvent,
        species::WashoutSpecies = WASHOUT_TRITIUM_IODINE,
)
    x ≥ 0 || throw(DomainError(x, "travel distance cannot be negative"))
    0 < u < Inf || throw(DomainError(u, "transport speed must be positive and finite"))
    f = precipitation_fraction(event)
    return 1 - f + f * exp(-washout_coefficients(event, species).low * x / u)
end

wet_depletion_factor(
    event::WashoutEvent,
    species::WashoutSpecies = WASHOUT_TRITIUM_IODINE,
) = exp(-washout_coefficients(event, species).low * event.duration)

# Airborne fraction after a travel time τx when a precipitation of coefficient Λ and
# duration t_w begins with a release of duration T_R, fixed over the path: the parcel
# released at τ is scavenged for min(τx, t_w − τ), the time until it arrives or the
# precipitation stops. With s = t_w − τ the release spans s ∈ [t_w − T_R, t_w]; the
# parcels with s < 0 leave after the precipitation and are not scavenged. For
# t_w ≤ T_R this is the most a precipitation lying within the release can leave
# airborne, against the 1 − f + f exp(−Λ τx) of wet_depletion_factor; for one longer
# than the release the most is that of t_w = T_R, a precipitation ending with it.
function _airborne_after_early_precipitation(τx::Real, Λ::Real, t_w::Real, T_R::Real)
    a, b = max(0.0, t_w - T_R), float(t_w)
    c = min(b, τx)          # exposure grows with s up to τx, then stays there
    rising = a < c ? (iszero(Λ) ? c - a : exp(-Λ * a) * -expm1(-Λ * (c - a)) / Λ) : 0.0
    flat = max(0.0, b - max(a, τx)) * exp(-Λ * τx)
    return (rising + flat + max(0.0, T_R - t_w)) / T_R
end

"""
    depletion_factor(x, site, class, nuclide; washout = nothing)

Total surviving airborne fraction over `x` metres: the product of the decay,
dry-deposition and washout factors.

They multiply because they act in sequence on the same material. With no
washout, no deposition and no decay each factor is one and so is the product,
which is the limit that fixes the composition. The 2021 code added the wet and
dry factors and so returned two in that limit.

`washout` is a [`WashoutEvent`](@ref); with the default of `nothing` the plume
travels in dry air.
"""
function depletion_factor(
        x::Real,
        site::AbstractSite,
        class::PasquillClass,
        nuclide::Nuclide;
        washout::Union{Nothing,WashoutEvent} = nothing,
)
    u = transport_wind_speed(site, class)
    f = decay_factor(x, u, nuclide) * dry_depletion_factor(x, site, class, nuclide)
    washout === nothing ||
        (f *= wet_depletion_factor(x, u, washout, nuclide.washout_species))
    return f
end
