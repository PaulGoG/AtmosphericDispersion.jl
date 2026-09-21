#=
Plume rise.

A plume leaves the stack with both momentum and buoyancy and keeps rising until
either ambient turbulence breaks it up or a stable stratification arrests it.
Each mechanism has a transitional law valid near the stack and a final rise
reached far from it; the model takes whichever mechanism dominates, or a
combined law when neither does.

Every correlation here is quoted for a *buoyant* plume in air whose stability
parameter is known. Where the 2021 code would have raised a domain error on a
negative base — an unstable atmosphere makes S negative, and S appears under
fractional powers — the inapplicable branches are dropped instead, which is
what the physics says: with no stable stratification there is no buoyancy
ceiling for it to impose.
=#

"""
    RiseCoefficients

The three plume-rise constants that published schemes disagree on, so that the
choice is made in the configuration rather than buried in the source.

  - `combined_buoyancy` — the denominator of the buoyancy term in the combined
    law. Briggs writes `3F x²/(2β²u³)` with an entrainment parameter `β = 0.6`,
    giving **0.72**; NSR-23 has **0.5**, which does not reduce to the two-thirds
    law as the momentum flux vanishes.
  - `neutral_momentum` — the coefficient of the neutral final momentum rise,
    `c w₀D/u`. Briggs (1969) Eq. 5.2, as implemented by EPA ISC3 Eq. (1-16),
    gives **3**. CNCAN NSR-23 has **1.5**, which is the momentum term of
    Holland's formula, `Δh = (w₀D/u)[1.5 + 2.68×10⁻³ p (T_s − T_a)D/T_s]` with
    `p` in mbar (Holland 1953; Turner, *Workbook of Atmospheric Dispersion
    Estimates*, 1970, Eq. 4.1), taken without its buoyancy term.
  - `stable_final` — the coefficient of the stable final buoyant rise,
    `c [F/(uS)]^(1/3)`. Briggs and the Handbook on Atmospheric Diffusion give
    **2.6**; NRC XOQDOQ writes **2.4**.

See [`BRIGGS_RISE`](@ref), [`XOQDOQ_RISE`](@ref) and [`NSR23_RISE`](@ref).
"""
struct RiseCoefficients
    combined_buoyancy::Float64
    neutral_momentum::Float64
    stable_final::Float64

    function RiseCoefficients(;
            combined_buoyancy::Real = 0.72,
            neutral_momentum::Real = 3.0,
            stable_final::Real = 2.6,
    )
        combined_buoyancy > 0 || throw(
            ArgumentError(
            "the combined-law buoyancy denominator must be positive, got $combined_buoyancy",
        ),
        )
        neutral_momentum ≥ 0 || throw(
            ArgumentError(
            "the neutral momentum coefficient cannot be negative, got $neutral_momentum",
        ),
        )
        stable_final ≥ 0 || throw(
            ArgumentError(
            "the stable rise coefficient cannot be negative, got $stable_final",
        ),
        )
        return new(combined_buoyancy, neutral_momentum, stable_final)
    end
end

"""
    BRIGGS_RISE

Briggs throughout, and the default: `2β² = 0.72`, `3 w₀D/u`, `2.6[F/(uS)]^(1/3)`.
"""
const BRIGGS_RISE = RiseCoefficients()

"""
    XOQDOQ_RISE

Briggs, with the stable coefficient NRC XOQDOQ (NUREG/CR-2919) writes as 2.4.
"""
const XOQDOQ_RISE = RiseCoefficients(stable_final = 2.4)

"""
    NSR23_RISE

The constants of CNCAN NSR-23, the Romanian normative, which the 2021 thesis
code followed: 0.5 in the combined law and Holland's 1.5 for the momentum rise.
`combined_buoyancy` here does not reduce to the two-thirds law.
"""
const NSR23_RISE = RiseCoefficients(combined_buoyancy = 0.5, neutral_momentum = 1.5)

"""
    buoyancy_transition_distance(F)

Distance `x₀` in metres over which a buoyant plume approaches its final rise,

    x₀ = 14 F^(5/8)   for F < 55 m⁴/s³
    x₀ = 34 F^(2/5)   otherwise

Final rise is attained by about `3.5 x₀`.
"""
function buoyancy_transition_distance(F::Real)
    F ≥ 0 || throw(DomainError(F, "the transition distance is defined for a buoyant plume"))
    return F < BUOYANCY_FLUX_BREAKPOINT ? 14 * F^(5 / 8) : 34 * F^(2 / 5)
end

"""
    BUOYANCY_FLUX_BREAKPOINT

Buoyancy flux in m⁴/s³ separating the two branches of the transition-distance
correlation.
"""
const BUOYANCY_FLUX_BREAKPOINT = 55.0

"""
    final_buoyant_rise(F, u, S, rise = BRIGGS_RISE)

Final rise in metres of a buoyant plume, from the buoyancy flux `F` in m⁴/s³,
the wind speed `u` in m/s at the release height, and the stability parameter
`S` in s⁻².

The smallest of the applicable limits is taken: the neutral limit, which
ambient turbulence sets and which always applies, and — only where the air is
stably stratified, `S > 0` — the stable and calm limits that stratification
imposes.
"""
function final_buoyant_rise(F::Real, u::Real, S::Real, rise::RiseCoefficients = BRIGGS_RISE)
    F ≥ 0 || throw(DomainError(F, "the rise correlations describe a buoyant plume"))
    u > 0 || throw(DomainError(u, "wind speed at the release height must be positive"))
    x₀ = buoyancy_transition_distance(F)
    neutral = 1.6 * F^(1 / 3) * (3.5 * x₀)^(2 / 3) / u
    S > 0 || return neutral
    stable = rise.stable_final * (F / (u * S))^(1 / 3)
    calm = 5.0 * F^(1 / 4) * S^(-3 / 8)
    return min(neutral, stable, calm)
end

"""
    buoyant_rise(x, F, u, S, rise = BRIGGS_RISE)

Rise in metres of a buoyant plume at downwind distance `x` metres: the
transitional law `1.6 F^(1/3) x^(2/3) / u` while the plume is still rising, and
the final rise beyond.
"""
function buoyant_rise(
        x::Real,
        F::Real,
        u::Real,
        S::Real,
        rise::RiseCoefficients = BRIGGS_RISE,
)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    final = final_buoyant_rise(F, u, S, rise)
    transitional = 1.6 * F^(1 / 3) * x^(2 / 3) / u
    x < 3.5 * buoyancy_transition_distance(F) && transitional ≤ final && return transitional
    return final
end

"""
    final_momentum_rise(Fₘ, w₀, D, u, S, rise = BRIGGS_RISE)

Final rise in metres driven by the exit momentum, from the momentum flux `Fₘ`
in m⁴/s², the exit velocity `w₀` in m/s, the stack diameter `D` in m, the wind
speed `u` in m/s and the stability parameter `S` in s⁻².

As for buoyancy, the stratification-limited branches enter only where `S > 0`.
"""
function final_momentum_rise(
        Fₘ::Real,
        w₀::Real,
        D::Real,
        u::Real,
        S::Real,
        rise::RiseCoefficients = BRIGGS_RISE,
)
    Fₘ ≥ 0 || throw(DomainError(Fₘ, "momentum flux cannot be negative"))
    u > 0 || throw(DomainError(u, "wind speed at the release height must be positive"))
    neutral = rise.neutral_momentum * w₀ * D / u
    S > 0 || return neutral
    calm = 4 * (Fₘ / S)^(1 / 4)
    stable = 1.5 * (Fₘ / u)^(1 / 3) * S^(-1 / 6)
    return min(neutral, calm, stable)
end

"""
    momentum_rise(x, Fₘ, w₀, D, u, S, rise = BRIGGS_RISE)

Rise in metres driven by exit momentum at downwind distance `x` metres, taking
the transitional law `1.89 (w₀² D / (u(w₀ + 3u)))^(2/3) x^(1/3)` until it
reaches the final rise.
"""
function momentum_rise(
        x::Real,
        Fₘ::Real,
        w₀::Real,
        D::Real,
        u::Real,
        S::Real,
        rise::RiseCoefficients = BRIGGS_RISE,
)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    final = final_momentum_rise(Fₘ, w₀, D, u, S, rise)
    transitional = 1.89 * (w₀^2 * D / (u * (w₀ + 3u)))^(2 / 3) * x^(1 / 3)
    return min(transitional, final)
end

"""
    combined_rise(x, F, Fₘ, w₀, u, S, D, rise = BRIGGS_RISE)

Rise in metres from the semi-empirical law combining momentum and buoyancy,

    Δh = 3^(1/3) [ Fₘ x / ((1/3 + u/w₀)² u²) + F x² / (c u³) ]^(1/3)

capped at the sum of the two final rises. Used where neither mechanism
dominates the other.

`c` is `rise.combined_buoyancy`. With Briggs' `c = 2β² = 0.72` the expression
reduces, as `Fₘ → 0`, to the two-thirds law `1.6 F^(1/3) x^(2/3) / u` that
[`buoyant_rise`](@ref) uses; with the 2021 value 0.5 it overshoots it by 13.6 %.
"""
function combined_rise(
        x::Real,
        F::Real,
        Fₘ::Real,
        w₀::Real,
        u::Real,
        S::Real,
        D::Real,
        rise::RiseCoefficients = BRIGGS_RISE,
)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    u > 0 || throw(DomainError(u, "wind speed at the release height must be positive"))
    w₀ > 0 || throw(DomainError(w₀, "exit velocity must be positive for the combined law"))
    final = final_momentum_rise(Fₘ, w₀, D, u, S, rise) + final_buoyant_rise(F, u, S, rise)
    transitional = 3^(1 / 3) *
                   (Fₘ * x / ((1 / 3 + u / w₀)^2 * u^2) +
                    F * x^2 / (rise.combined_buoyancy * u^3))^(
        1 / 3
    )
    return min(transitional, final)
end

"""
    MECHANISM_BALANCE_TOLERANCE

Relative difference between the two final rises below which neither mechanism
is taken to dominate and the combined law is used instead.
"""
const MECHANISM_BALANCE_TOLERANCE = 0.1

# The mechanism selection, shared by both public methods of plume_rise.
function _plume_rise(
        x::Real,
        F::Real,
        Fₘ::Real,
        S::Real,
        w₀::Real,
        D::Real,
        u::Real,
        rise::RiseCoefficients,
)
    buoyant = final_buoyant_rise(F, u, S, rise)
    momentum = final_momentum_rise(Fₘ, w₀, D, u, S, rise)
    total = buoyant + momentum
    balanced = iszero(total) ||
               2 * abs(buoyant - momentum) / total ≤ MECHANISM_BALANCE_TOLERANCE
    balanced && return combined_rise(x, F, Fₘ, w₀, u, S, D, rise)
    momentum > buoyant && return momentum_rise(x, Fₘ, w₀, D, u, S, rise)
    return buoyant_rise(x, F, u, S, rise)
end

"""
    plume_rise(x, source, atmosphere, u, rise = BRIGGS_RISE)

Rise of the plume above the release height, in metres, at downwind distance `x`
metres, with `u` the wind speed at the release height in m/s.

The dominant mechanism is chosen by comparing the two final rises: where they
agree to within [`MECHANISM_BALANCE_TOLERANCE`](@ref) in relative terms neither
dominates and the combined law applies, otherwise the larger mechanism carries
the plume and its law is used alone. `rise` selects the
[`RiseCoefficients`](@ref).
"""
function plume_rise(
        x::Real,
        source::StackSource,
        atmosphere::Atmosphere,
        u::Real,
        rise::RiseCoefficients = BRIGGS_RISE,
)
    F = buoyancy_flux(source, atmosphere)
    Fₘ = momentum_flux(source, atmosphere)
    S = stability_parameter(atmosphere)
    return _plume_rise(x, F, Fₘ, S, source.exit_velocity, source.diameter, u, rise)
end
