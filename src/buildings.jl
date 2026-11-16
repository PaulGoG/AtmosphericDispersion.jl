"""
    Building(; east, north, height, frontal_area)

A structure near the stack, positioned relative to it in the geographic frame.

- `east`, `north` — displacement from the stack, m
- `height` — building height above ground, m
- `frontal_area` — cross-sectional area presented to the wind, m²
"""
struct Building
    east::Float64
    north::Float64
    height::Float64
    frontal_area::Float64

    function Building(; east::Real, north::Real, height::Real, frontal_area::Real)
        height ≥ 0 ||
            throw(ArgumentError("building height cannot be negative, got $height m"))
        frontal_area ≥ 0 ||
            throw(ArgumentError("frontal area cannot be negative, got $frontal_area m²"))
        return new(east, north, height, frontal_area)
    end
end

"""
    distance(building)

Horizontal distance of `building` from the stack, in metres.
"""
distance(b::Building) = hypot(b.east, b.north)

"""
    WAKE_INFLUENCE_RADII

Multiple of a building's own height within which it is taken to influence the
plume; buildings further from the stack are ignored.
"""
const WAKE_INFLUENCE_RADII = 3.0

"""
    DEFAULT_WAKE_COEFFICIENT

Coefficient `C` of the wake broadening of the dispersion parameters,
`√(σ² + C A/π)`. Setting it to zero disables the wake broadening; the release
height is still lowered by [`wake_height`](@ref).

The default is one, which IAEA Safety Reports Series No. 19 [IAEA2001](@cite)
Eq. (6) writes as `Σ_z = (σ_z² + A_B/π)^(1/2)`, for the vertical parameter
only, and the German AVV zu §47 StrlSchV [AVV2012](@cite) Eqs. (4.31)/(4.32)
as `√(σ² + I_G²/π)`. This package broadens both parameters, as CNCAN NSR-23
[CNCAN2004](@cite) Eq. (9) does. NSR-23 gives `C` between 0.5 and 2.0 and
recommends the value within that interval which gives the largest dilution
factor; the 2021 thesis code used 1.5. NRC Regulatory Guide 1.111
[NRC1977](@cite) Eq. (9) is a third convention again, applying 0.5 to the building
height rather than its area.

See [`NSR23_WAKE_COEFFICIENT`](@ref) to restore the 2021 value.
"""
const DEFAULT_WAKE_COEFFICIENT = 1.0

"""
    NSR23_WAKE_COEFFICIENT

The wake coefficient of the 2021 thesis code, kept so that its results can be
reproduced. It lies inside the interval 0.5–2.0 that CNCAN NSR-23 gives; see
[`DEFAULT_WAKE_COEFFICIENT`](@ref).
"""
const NSR23_WAKE_COEFFICIENT = 1.5

"""
    MAXIMUM_WAKE_BROADENING

Largest factor by which the building wake may widen a dispersion parameter, √3.
NRC XOQDOQ [Sagendorf1982](@cite), after Regulatory Guide 1.111, bounds the
"maximum allowable building-wake dilution" by `3σ_yσ_z` in its point form,
Eq. (8), and by `√3 σ_z` in its sector form, Eq. (11); CNCAN NSR-23
[CNCAN2004](@cite) §5 recommends that the reduction of concentration by the
wake "is not higher than 3". Limiting each parameter to √3 reproduces the
sector form exactly and keeps the reduction of a centreline concentration
within 3. Where the two parameters widen unequally, as they do when
`σ_y ≫ σ_z`, it is stricter than the product bound of Eq. (8), which errs
towards the higher concentration.
"""
const MAXIMUM_WAKE_BROADENING = sqrt(3.0)

"""
    BuildingEnvelope(buildings = Building[]; wake_coefficient = DEFAULT_WAKE_COEFFICIENT,
                     broadening_limit = MAXIMUM_WAKE_BROADENING)

The building that governs the wake of a stack, reduced once to a height and a
frontal area, and the rules of the wake correction: the coefficient `C` of
[`wake_broadened`](@ref) and the largest factor `broadening_limit` by which it
may widen a dispersion parameter (`Inf` for no limit, as in IAEA SRS-19).

Buildings within [`WAKE_INFLUENCE_RADII`](@ref) of their own height from the
stack are in its vicinity; the rest are ignored. Of those, the tallest governs,
and of equally tall ones the one with the larger frontal area. IAEA SRS-19
[IAEA2001](@cite) §3.3 takes "the building that most influences the resulting
plume dispersion", noting that a much larger building nearby exerts more
influence than the one the release comes from; CNCAN NSR-23 [CNCAN2004](@cite)
§5 takes the nearest building in the vicinity. The 2021 code averaged the
heights and areas with weights inversely proportional to distance, so that a
small structure near the stack could cancel the wake of a large one.

The reduction is performed at construction. The 2021 code recomputed it inside
the dispersion-parameter correction, which is called once per grid point per
stability class, so the building list, which does not change during a run, was
rescanned at each evaluation.

An empty envelope has zero height and area, which makes every building
correction downstream the identity.
"""
struct BuildingEnvelope
    height::Float64
    frontal_area::Float64
    wake_coefficient::Float64
    broadening_limit::Float64

    function BuildingEnvelope(
            buildings::AbstractVector{Building} = Building[];
            wake_coefficient::Real = DEFAULT_WAKE_COEFFICIENT,
            broadening_limit::Real = MAXIMUM_WAKE_BROADENING,
    )
        wake_coefficient ≥ 0 || throw(
            ArgumentError("the wake coefficient cannot be negative, got $wake_coefficient"),
        )
        broadening_limit ≥ 1 || throw(
            ArgumentError("the broadening limit cannot be below one, got $broadening_limit"),
        )
        h, a = 0.0, 0.0
        for b in buildings
            d = distance(b)
            d > 0 || throw(
                ArgumentError("a building cannot sit at the foot of the stack itself"),
            )
            d ≤ WAKE_INFLUENCE_RADII * b.height || continue
            (b.height, b.frontal_area) > (h, a) && ((h, a) = (b.height, b.frontal_area))
        end
        return new(h, a, float(wake_coefficient), float(broadening_limit))
    end
end

"""
    equivalent_height(envelope)

Height in metres of the building that governs the wake. Zero when no building is
close enough to influence the plume.
"""
equivalent_height(envelope::BuildingEnvelope) = envelope.height

"""
    equivalent_area(envelope)

Frontal area in m² of the building that governs the wake.
"""
equivalent_area(envelope::BuildingEnvelope) = envelope.frontal_area

"""
    wake_broadened(σ, H, envelope)

Dispersion parameter in metres, broadened by the building wake.

A plume released well above the buildings (higher than
`2.5 × equivalent_height`) is unaffected. One released below their tops is
fully mixed into the wake and takes the broadened value
`min(√(σ² + C A / π), L σ)`, with `L` the `broadening_limit` of the envelope.
Between the two the correction is interpolated linearly in release height, as
in CNCAN NSR-23 [CNCAN2004](@cite) Eq. (10).

With an empty envelope this is the identity, since the height and area are
then zero.
"""
function wake_broadened(σ::Real, H::Real, envelope::BuildingEnvelope)
    h = equivalent_height(envelope)
    H ≥ 2.5 * h && return float(σ)
    broadened = min(
        sqrt(σ^2 + envelope.wake_coefficient * equivalent_area(envelope) / π),
        envelope.broadening_limit * σ,
    )
    H < h && return broadened
    return broadened - (H - h) / (1.5 * h) * (broadened - σ)
end

"""
    BuildingZone

Where a receptor lies relative to the building that most disturbs the flow,
after IAEA Safety Reports Series No. 19 [IAEA2001](@cite) §3.3, which bounds
the three zones by the release height `H`, the building height `H_B`, its
frontal area `A_B` and the downwind distance `x`:

  - `DISPLACEMENT_ZONE` — `H > 2.5 H_B`: the flow is deflected around the
    building and dispersion is undisturbed.
  - `WAKE_ZONE` — `H ≤ 2.5 H_B` and `x > 2.5 √A_B`: the disturbed wake downwind
    of the building, where the plume is broadened by [`wake_broadened`](@ref).
  - `CAVITY_ZONE` — `H ≤ 2.5 H_B` and `x ≤ 2.5 √A_B`: the recirculating cavity
    in the lee of the building, where the Gaussian plume does not apply and
    [`dilution_cavity`](@ref) or [`dilution_cavity_wall`](@ref) is used.
"""
@enum BuildingZone::UInt8 begin
    DISPLACEMENT_ZONE = 1
    WAKE_ZONE = 2
    CAVITY_ZONE = 3
end

"""
    CAVITY_EXTENT_FACTOR

Multiple of `√A_B`, the square root of the building's frontal area, that the
cavity zone extends downwind of it: 2.5, IAEA SRS-19 §3.3.
"""
const CAVITY_EXTENT_FACTOR = 2.5

"""
    building_zone(H, x, envelope)

The [`BuildingZone`](@ref) of a receptor at downwind distance `x` metres for a
plume at effective height `H` metres over the [`BuildingEnvelope`](@ref). An
empty envelope places every receptor in the displacement zone.
"""
function building_zone(H::Real, x::Real, envelope::BuildingEnvelope)
    H ≥ 0 || throw(DomainError(H, "effective height cannot be below ground"))
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    h = equivalent_height(envelope)
    (iszero(h) || H > 2.5 * h) && return DISPLACEMENT_ZONE
    x > CAVITY_EXTENT_FACTOR * sqrt(equivalent_area(envelope)) && return WAKE_ZONE
    return CAVITY_ZONE
end
