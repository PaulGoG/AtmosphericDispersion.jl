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
plume. A building further from the stack than this plays no part.
"""
const WAKE_INFLUENCE_RADII = 3.0

"""
    DEFAULT_WAKE_COEFFICIENT

Coefficient `C` of the wake broadening of the dispersion parameters,
`√(σ² + C A/π)`. Setting it to zero disables the building correction entirely.

**One**, which is what IAEA Safety Reports Series No. 19 [IAEA2001](@cite) Eq. (6) writes
as `Σ_z = (σ_z² + A_B/π)^(1/2)` and the German AVV zu §47 StrlSchV [AVV2012](@cite)
Eqs. (4.31)/(4.32) as `√(σ² + I_G²/π)`. The 2021 thesis code used 1.5, following
the Romanian normative it was written against, and no source outside that
normative was found for it; NRC Regulatory Guide 1.111 [NRC1977](@cite) Eq. (9) is a third
convention again, applying 0.5 to the building *height* rather than its area.

See [`NSR23_WAKE_COEFFICIENT`](@ref) to restore the 2021 value.
"""
const DEFAULT_WAKE_COEFFICIENT = 1.0

"""
    NSR23_WAKE_COEFFICIENT

The wake coefficient of the Romanian normative the 2021 thesis followed, kept so
its results can be reproduced. See [`DEFAULT_WAKE_COEFFICIENT`](@ref).
"""
const NSR23_WAKE_COEFFICIENT = 1.5

"""
    BuildingEnvelope(buildings = Building[]; wake_coefficient = DEFAULT_WAKE_COEFFICIENT)

The collective effect of the buildings around a stack, reduced once to a single
equivalent height and frontal area.

Buildings within [`WAKE_INFLUENCE_RADII`](@ref) of their own height from the
stack contribute; the rest are ignored. Contributions are averaged with weights
inversely proportional to distance, so nearer structures dominate.

The reduction is performed at construction. The 2021 code recomputed it inside
the dispersion-parameter correction, which is called once per grid point per
stability class, so the whole building list was rescanned for every evaluation
of a field that never changes.

An empty envelope has zero equivalent height and area, which makes every
building correction downstream the identity.
"""
struct BuildingEnvelope
    height::Float64
    frontal_area::Float64
    wake_coefficient::Float64

    function BuildingEnvelope(
            buildings::AbstractVector{Building} = Building[];
            wake_coefficient::Real = DEFAULT_WAKE_COEFFICIENT,
    )
        wake_coefficient ≥ 0 || throw(
            ArgumentError("the wake coefficient cannot be negative, got $wake_coefficient"),
        )
        weight = 0.0
        h = 0.0
        a = 0.0
        for b in buildings
            d = distance(b)
            d > 0 || throw(
                ArgumentError("a building cannot sit at the foot of the stack itself"),
            )
            d ≤ WAKE_INFLUENCE_RADII * b.height || continue
            w = inv(d)
            weight += w
            h += w * b.height
            a += w * b.frontal_area
        end
        iszero(weight) && return new(0.0, 0.0, float(wake_coefficient))
        return new(h / weight, a / weight, float(wake_coefficient))
    end
end

"""
    equivalent_height(envelope)

Weighted equivalent building height, in metres. Zero when no building is close
enough to influence the plume.
"""
equivalent_height(envelope::BuildingEnvelope) = envelope.height

"""
    equivalent_area(envelope)

Weighted equivalent frontal area, in m².
"""
equivalent_area(envelope::BuildingEnvelope) = envelope.frontal_area

"""
    wake_broadened(σ, H, envelope)

Dispersion parameter in metres, broadened by the building wake.

A plume released well above the buildings — higher than
`2.5 × equivalent_height` — is unaffected. One released below their tops is
fully mixed into the wake and takes the broadened value
`√(σ² + C A / π)`. Between the two the correction is interpolated linearly in
release height.

With an empty envelope this is the identity, since every release height clears
a ceiling of zero.
"""
function wake_broadened(σ::Real, H::Real, envelope::BuildingEnvelope)
    h = equivalent_height(envelope)
    H ≥ 2.5 * h && return float(σ)
    broadened = sqrt(σ^2 + envelope.wake_coefficient * equivalent_area(envelope) / π)
    H < h && return broadened
    return broadened - (H - h) / (1.5 * h) * (broadened - σ)
end
