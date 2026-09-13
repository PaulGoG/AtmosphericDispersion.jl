"""
    DepositionVelocity(low, high)

Bracketing dry deposition velocities in m/s for a nuclide over a given surface.

The two are used in opposite places, and deliberately so. Depletion of the
airborne plume takes `low`, because removing less material leaves more in the
air and overestimates the inhalation pathway. Ground deposition takes `high`,
because depositing more material overestimates the ground-shine and ingestion
pathways. Each choice is conservative for the endpoint it feeds.
"""
struct DepositionVelocity
    low::Float64
    high::Float64

    function DepositionVelocity(low::Real, high::Real)
        0 ≤ low ≤ high || throw(
            ArgumentError(
                "deposition velocities must satisfy 0 ≤ low ≤ high, got low = $low, high = $high m/s",
            ),
        )
        return new(low, high)
    end
end

"""
    Nuclide(; name, decay_constant, deposition_velocity)

A released species: what it is called, how fast it decays, and how readily it
deposits.

- `decay_constant` — radioactive decay constant `λ`, s⁻¹
- `deposition_velocity` — [`DepositionVelocity`](@ref) over the receptor surface

The dispersion core is otherwise independent of the species; a nuclide enters
only through depletion and deposition.
"""
struct Nuclide
    name::String
    decay_constant::Float64
    deposition_velocity::DepositionVelocity

    function Nuclide(;
        name::AbstractString,
        decay_constant::Real,
        deposition_velocity::DepositionVelocity,
    )
        decay_constant ≥ 0 || throw(
            ArgumentError("decay constant cannot be negative, got $decay_constant s⁻¹"),
        )
        return new(String(name), float(decay_constant), deposition_velocity)
    end
end

Base.show(io::IO, n::Nuclide) = print(io, "Nuclide(", n.name, ")")

"""
    TRITIUM_DECAY_CONSTANT

Decay constant of tritium, 1.784e-9 s⁻¹, a half-life of about 12.3 years.
"""
const TRITIUM_DECAY_CONSTANT = 1.784e-9

"""
    TRITIATED_WATER

Tritium as HTO, the chemical form a heavy-water reactor releases in quantity
and the one that deposits readily.
"""
const TRITIATED_WATER = Nuclide(;
    name = "HTO",
    decay_constant = TRITIUM_DECAY_CONSTANT,
    deposition_velocity = DepositionVelocity(0.4e-2, 0.8e-2),
)

"""
    TRITIUM_GAS

Tritium as HT, an order of magnitude less depositing than [`TRITIATED_WATER`](@ref).
"""
const TRITIUM_GAS = Nuclide(;
    name = "HT",
    decay_constant = TRITIUM_DECAY_CONSTANT,
    deposition_velocity = DepositionVelocity(0.04e-2, 0.05e-2),
)
