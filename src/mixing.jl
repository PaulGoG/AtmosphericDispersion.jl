#=
The mixing layer.

A plume does not disperse upwards without limit. Turbulent mixing is confined
below a capping inversion, and once the vertical dispersion parameter approaches
that depth the plume is trapped between the ground and the lid, reflecting off
both. Beyond that the vertical profile becomes uniform and the concentration
falls as 1/x rather than 1/(x σ_z).

Two regulatory conventions exist for it, and they differ in more than the
depths. HPA-RPD-058 §3.2.2.1, after NRPB-R91, caps every stability category at
a tabulated depth. EPA ISC3 (User's Guide vol. II, §1.1.6.1) takes measured
mixing heights, applies them to the unstable and neutral categories only, and
lets a plume whose effective height exceeds the mixing height leave the layer
altogether. `MixingLayer` carries what distinguishes them, so either can be
stated rather than approximated.
=#

"""
    LidRule

What happens to a plume whose effective height lies above the top of the mixing
layer.

  - `RISE_INHIBITED` — the capping inversion arrests the rise, and the plume is
    dispersed from the top of the layer: the effective height is taken as
    `min(H, A)`. NRPB-R157 §B2.3: "plume rise will be inhibited by a capping
    inversion to the mixing layer. If a plume rises into such an inversion the
    amount of material in the mixing layer, and hence ground-level
    concentration, will be reduced." Keeping all of it in the layer is therefore
    the conservative reading, it is continuous in `H`, and it is what
    HPA-RPD-058 Table 3.7 tabulates for a 100 m release under a 100 m lid. The
    default. A source that stands physically above the lid is treated the same
    way.
  - `FULL_PENETRATION` — EPA ISC3: "if the effective stack height exceeds the
    mixing height, the plume is assumed to fully penetrate the elevated
    inversion and the ground-level concentration is set equal to zero." Nothing
    reaches the layer, so nothing deposits from it either.
"""
@enum LidRule begin
    RISE_INHIBITED = 1
    FULL_PENETRATION = 2
end

"""
    MixingLayer(depths; above_lid = RISE_INHIBITED)
    MixingLayer(depth; above_lid = RISE_INHIBITED)

Depth of the mixing layer in metres for each Pasquill class A to F, and the
[`LidRule`](@ref) applied to a plume above it. `Inf` leaves a class unbounded. A
single `depth` applies to every class.

The depth is a site measurement. Where it has not been measured,
[`MIXING_TABULATED`](@ref) carries typical values; where it has, pass them:

```julia
MixingLayer((1500.0, 1200.0, 1000.0, 900.0, 300.0, 120.0))
```

The ISC3 convention, which treats stable air as unbounded and lets a plume
through the lid, is

```julia
MixingLayer((1500.0, 1200.0, 1000.0, 900.0, Inf, Inf); above_lid = FULL_PENETRATION)
```
"""
struct MixingLayer
    depths::NTuple{6,Float64}
    above_lid::LidRule

    function MixingLayer(depths; above_lid::LidRule = RISE_INHIBITED)
        length(depths) == 6 || throw(
            ArgumentError(
                "one mixing depth per Pasquill class A to F is required, got $(length(depths))",
            ),
        )
        all(d -> d isa Real && d > 0, depths) || throw(
            ArgumentError("mixing depths must be positive, `Inf` for none, got $depths"),
        )
        return new(ntuple(i -> Float64(depths[i]), 6), above_lid)
    end
end

MixingLayer(depth::Real; above_lid::LidRule = RISE_INHIBITED) =
    MixingLayer(ntuple(_ -> depth, 6); above_lid)

"""
    MIXING_TABULATED

HPA-RPD-058 Table 3.5(a), "Typical values of wind speed and depth of mixing
layer for use when measured values are not available", attributed there to
Clarke (1979) and Jones (1980): A 1300, B 900, C 850, D 800, E 400, F 100 m. The
default. The table's category G is not carried, this package having none.
"""
const MIXING_TABULATED = MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0, 100.0))

"""
    MIXING_UNBOUNDED

No lid in any class: the plume disperses upwards without limit.
"""
const MIXING_UNBOUNDED = MixingLayer(Inf)

"""
    RECOMMENDED_MIXING_DEPTH

The single depth in metres HPA-RPD-058 §3.2.2.2.3 recommends "for all
conditions" where the real one is unknown; `MixingLayer(RECOMMENDED_MIXING_DEPTH)`.
"""
const RECOMMENDED_MIXING_DEPTH = 800.0

"""
    mixing_depth(class, layer = MIXING_TABULATED)

Depth of the mixing layer in metres for the given Pasquill `class`, `Inf` where
the class is unbounded.
"""
mixing_depth(class::PasquillClass, layer::MixingLayer = MIXING_TABULATED) =
    layer.depths[classindex(class)]

# Ratio Σ_z/A at which the evaluation changes representation; see `_lidded_profile`.
const _REPRESENTATION_CROSSOVER = 0.7

# The vertical profile between two reflecting planes at 0 and A, for 0 ≤ z, H ≤ A,
# normalised to unit integral over [0, A].
#
# HPA-RPD-058 Eq. (3.4) writes it as a sum over image sources at 2sA ± H and
# notes that the series "can be summed to any prescribed accuracy". It converges
# fast while Σ_z is small against A and slowly once it is not, where the same
# function is better written as its cosine series,
#
#     (1/A) [1 + 2 Σ_k exp(−k²π²Σ_z²/2A²) cos(kπz/A) cos(kπH/A)]
#
# which converges fast exactly there and whose leading term is the uniform
# profile 1/A of Eq. (3.5). Each is used on its own side of the crossover, to
# double precision: the neglected image terms are below exp(−50) and the
# neglected harmonics below exp(−38).
function _lidded_profile(z::Float64, H::Float64, Σz::Float64, A::Float64)
    if Σz < _REPRESENTATION_CROSSOVER * A
        total = 0.0
        for s = -4:4
            total += exp(-(z - H - 2s * A)^2 / (2Σz^2)) + exp(-(z + H - 2s * A)^2 / (2Σz^2))
        end
        return total / (sqrt(2π) * Σz)
    end
    total = 1.0
    for k = 1:4
        total += 2 * exp(-(k * π * Σz / A)^2 / 2) * cos(k * π * z / A) * cos(k * π * H / A)
    end
    return total / A
end

"""
    vertical_factor(z, H, Σz, A, rule = RISE_INHIBITED)
    vertical_factor(z, H, Σz, layer, class)

The vertical part of the Gaussian plume at receptor height `z` for a release at
effective height `H` with vertical dispersion `Σz`, under a mixing layer of
depth `A` metres, in m⁻¹.

With `A = Inf` it is the ordinary ground-reflected Gaussian,

    [exp(−(z − H)²/2Σ_z²) + exp(−(z + H)²/2Σ_z²)] / (√(2π) Σ_z)

Under a lid it is the profile between two reflecting planes, HPA-RPD-058
Eq. (3.4) summed to convergence, which goes over to the uniform `1/A` of
Eq. (3.5) as `Σ_z` passes the depth. The activity between the ground and the
lid integrates to one at every `Σ_z`, and the factor is zero above the lid.

A plume above the lid, `H > A`, is treated by `rule`; see [`LidRule`](@ref).
"""
function vertical_factor(
    z::Real,
    H::Real,
    Σz::Real,
    A::Real,
    rule::LidRule = RISE_INHIBITED,
)
    Σz > 0 || throw(DomainError(Σz, "the vertical dispersion parameter must be positive"))
    A > 0 || throw(DomainError(A, "the mixing depth must be positive"))
    z ≥ 0 || throw(DomainError(z, "receptor height cannot be below ground"))
    H ≥ 0 || throw(DomainError(H, "effective height cannot be below ground"))
    if !isfinite(A)
        return (exp(-(z - H)^2 / (2Σz^2)) + exp(-(z + H)^2 / (2Σz^2))) / (sqrt(2π) * Σz)
    end
    z > A && return 0.0
    if H > A
        rule == FULL_PENETRATION && return 0.0
        H = A
    end
    return _lidded_profile(Float64(z), Float64(H), Float64(Σz), Float64(A))
end

vertical_factor(z::Real, H::Real, Σz::Real, layer::MixingLayer, class::PasquillClass) =
    vertical_factor(z, H, Σz, mixing_depth(class, layer), layer.above_lid)

"""
    crosswind_integrated_factor(H, Σz, A, rule = RISE_INHIBITED)
    crosswind_integrated_factor(H, Σz, layer, class)

The vertical part of the crosswind-integrated plume, in m⁻¹:
[`vertical_factor`](@ref) at ground level. Without a lid the two image terms
coincide and it is exactly the `√(2/π) exp(−H²/2Σ_z²)/Σ_z` of the
sector-averaged form; under a lid the plume has filled it tends to `1/A`.
"""
crosswind_integrated_factor(H::Real, Σz::Real, A::Real, rule::LidRule = RISE_INHIBITED) =
    vertical_factor(0.0, H, Σz, A, rule)

crosswind_integrated_factor(H::Real, Σz::Real, layer::MixingLayer, class::PasquillClass) =
    vertical_factor(0.0, H, Σz, layer, class)
