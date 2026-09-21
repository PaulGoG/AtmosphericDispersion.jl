#=
The mixing layer.

A plume does not disperse upwards without limit. Turbulent mixing is confined
below a capping inversion, and once the vertical dispersion parameter approaches
that depth the plume is trapped between the ground and the lid, reflecting off
both. Beyond that the vertical profile becomes uniform and the concentration
falls as 1/x rather than 1/(x σ_z).

Its absence is what limited agreement with published depletion tables beyond
about 20 km, where an unbounded plume keeps diluting vertically and a real one
cannot.
=#

"""
    MixingLayer

How the top of the mixing layer is treated.

  - `MIXING_TABULATED` — the depth is taken per stability class from
    [`mixing_depth`](@ref), the default.
  - `MIXING_UNBOUNDED` — no lid, the plume disperses upwards without limit. This
    is what the 2021 code did and what this package did before the lid existed.
  - `MIXING_UNIFORM_800` — a single 800 m depth for every class, which is what
    the EC group behind HPA-RPD-058 §3.2.2.2.3 recommend "for all conditions" as
    a compromise when the real depth is unknown.
"""
@enum MixingLayer begin
    MIXING_TABULATED = 1
    MIXING_UNBOUNDED = 2
    MIXING_UNIFORM_800 = 3
end

# HPA-RPD-058 Table 3.5(a), "Typical values of wind speed and depth of mixing
# layer for use when measured values are not available", Pasquill/Smith/Hosker
# scheme, attributed there to Clarke (1979) and Jones (1980). The G row, 100 m,
# is not carried because this package has no category G.
const _MIXING_DEPTH = (1300.0, 900.0, 850.0, 800.0, 400.0, 100.0)

"""
    RECOMMENDED_MIXING_DEPTH

The single depth in metres recommended for all conditions where the real one is
unknown, HPA-RPD-058 §3.2.2.2.3.
"""
const RECOMMENDED_MIXING_DEPTH = 800.0

"""
    mixing_depth(class, layer = MIXING_TABULATED)

Depth of the mixing layer in metres for the given Pasquill `class`, or `Inf`
where the plume is taken as unbounded.

The tabulated values fall from 1300 m in very unstable air to 100 m in very
stable air, which is the physical ordering: convection deepens the mixed layer
and stable stratification suppresses it.
"""
function mixing_depth(class::PasquillClass, layer::MixingLayer = MIXING_TABULATED)
    layer == MIXING_UNBOUNDED && return Inf
    layer == MIXING_UNIFORM_800 && return RECOMMENDED_MIXING_DEPTH
    return _MIXING_DEPTH[classindex(class)]
end

"""
    UNIFORM_MIXING_RATIO

Ratio of `σ_z` to the mixing depth beyond which the vertical profile is taken as
uniform rather than summed over images.

HPA-RPD-058 §3.2.2.1 puts the transition at "when the value of the vertical
dispersion coefficient becomes greater than the depth of the mixing layer", so
the switch is made at 1.0. The two descriptions do not quite meet there: with
the image sum truncated at `|s| = 1` the profile at `Σ_z = A` still departs from
`1/A`, by 1.4 % at the ground for a low release and by up to 4 % at the lid for
a release at the lid, so the factor steps by that much across the switch.
"""
const UNIFORM_MIXING_RATIO = 1.0

"""
    vertical_factor(z, H, Σz, A)

The vertical part of the Gaussian plume at receptor height `z` for a release at
effective height `H` with vertical dispersion `Σz`, under a mixing layer of
depth `A` metres, in m⁻¹.

Below the lid this is the sum over the ground and lid images of HPA-RPD-058
Eq. (3.4),

    Σ_s exp(−(z − H − 2sA)²/2Σ_z²) + exp(−(z + H − 2sA)²/2Σ_z²)

divided by `√(2π) Σ_z`, truncated at `|s| = 1` as the report prescribes. Once
`Σ_z` reaches the depth the profile is uniform and it becomes `1/A`, which is
Eq. (3.5).

Above the lid the factor is zero: the material is confined to `0 ≤ z ≤ A`.
Integrating over that interval returns the whole release to better than 5 parts
in 10⁴ while `Σ_z ≤ 0.6 A`. Beyond that the truncation of the image sum starts
to tell: just below the switch to the uniform profile, which is exact again,
the integral is 0.997 for a release at a tenth of the depth and 0.977 for one at
the lid itself.

**A release strictly above the lid is not treated.** HPA-RPD-058's Diagram 3.1
places the source below the inversion, and a plume released above one is
decoupled from the ground until the inversion breaks — fumigation, which is a
different model and not one this package implements. Where `H > A` the lid is
therefore ignored and the unbounded form used, which is what happens in class F
for the reference case: its tabulated depth is 100 m and plume rise puts the
release at 103 m. Trapping such a plume against a lid it is already above would
roughly double the ground-level concentration on no physical grounds.

A release sitting *exactly* at the lid is still capped. That boundary is not a
matter of taste: HPA-RPD-058's own Table 3.7 tabulates a 100 m release in
category F, whose depth is 100 m, and only the capped form reproduces it — 0.357
and 0.089 at 50 and 100 km against a published 0.33 and 0.083, where the
uncapped form gives 0.599 and 0.311.

With `A = Inf` the sum collapses to the two unbounded terms and the result is
the ordinary ground-reflected Gaussian.
"""
function vertical_factor(z::Real, H::Real, Σz::Real, A::Real)
    Σz > 0 || throw(DomainError(Σz, "the vertical dispersion parameter must be positive"))
    isfinite(A) && A ≤ 0 && throw(DomainError(A, "the mixing depth must be positive"))
    # The lid applies only to a release below it; see the note above.
    capped = isfinite(A) && H ≤ A
    if capped
        # Nothing crosses the lid. Without this the profile would carry material
        # above it and the released activity would not be conserved.
        z > A && return 0.0
        Σz ≥ UNIFORM_MIXING_RATIO * A && return 1 / A
    end
    total = exp(-(z - H)^2 / (2Σz^2)) + exp(-(z + H)^2 / (2Σz^2))
    if capped
        for s in (-1, 1)
            total += exp(-(z - H - 2s * A)^2 / (2Σz^2)) + exp(-(z + H - 2s * A)^2 / (2Σz^2))
        end
    end
    return total / (sqrt(2π) * Σz)
end

"""
    crosswind_integrated_factor(H, Σz, A)

The vertical part of the crosswind-integrated plume, in m⁻¹: [`vertical_factor`](@ref)
at ground level, where the two image terms coincide and their sum is exactly the
`√(2/π) exp(−H²/2Σ_z²)/Σ_z` of the sector-averaged form.

Reduces to `1/A` under a lid the plume has filled, which is Eq. (3.5).
"""
crosswind_integrated_factor(H::Real, Σz::Real, A::Real) = vertical_factor(0.0, H, Σz, A)
