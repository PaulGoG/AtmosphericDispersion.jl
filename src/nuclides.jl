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
    Nuclide(; name, decay_constant, deposition_velocity, washout_species = WASHOUT_TRITIUM_IODINE)

A released species: what it is called, how fast it decays, and how readily it
deposits.

- `decay_constant` — radioactive decay constant `λ`, s⁻¹
- `deposition_velocity` — [`DepositionVelocity`](@ref) over the receptor surface
- `washout_species` — the [`WashoutSpecies`](@ref) row of the washout table the nuclide falls under

The dispersion core is otherwise independent of the species; a nuclide enters
only through depletion and deposition.
"""
struct Nuclide
    name::String
    decay_constant::Float64
    deposition_velocity::DepositionVelocity
    washout_species::WashoutSpecies

    function Nuclide(;
            name::AbstractString,
            decay_constant::Real,
            deposition_velocity::DepositionVelocity,
            washout_species::WashoutSpecies = WASHOUT_TRITIUM_IODINE,
    )
        decay_constant ≥ 0 || throw(
            ArgumentError("decay constant cannot be negative, got $decay_constant s⁻¹"),
        )
        return new(
            String(name),
            float(decay_constant),
            deposition_velocity,
            washout_species,
        )
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

Deposition velocity 0.4–0.8 × 10⁻² m/s, from the notes to CNCAN NSR-23 Table 6,
which attribute it to experimental measurement — [Murphy1993](@citet) — and
attach the condition that the tropopause be taken at 12–15 km. The lower bound
agrees with the values in common use elsewhere: the MACCS2 default of 0.5,
0.42 measured at the Savannah River Site, and 0.392–0.444 from AECL.

**That HTO deposits at all is a divergence from other guidance**, and a
deliberate one. IAEA SRS-19 [IAEA2001](@cite) §3.9, EUR 15760 [Simmonds1995](@cite)
§3.2 and HPA-RPD-058 [SmithSimmonds2009](@cite) §3.2.2.3 each assign tritium a deposition velocity of **zero** and handle it by specific
activity instead. This package follows the normative the work was done under;
set the velocities to zero to follow the others.
"""
const TRITIATED_WATER = Nuclide(;
    name = "HTO",
    decay_constant = TRITIUM_DECAY_CONSTANT,
    deposition_velocity = DepositionVelocity(0.4e-2, 0.8e-2),
)

"""
    TRITIUM_GAS

Tritium as HT, an order of magnitude less depositing than [`TRITIATED_WATER`](@ref).

Deposition velocity 0.04–0.05 × 10⁻² m/s, the same source, which puts HT an
order of magnitude below HTO under the same conditions.
"""
const TRITIUM_GAS = Nuclide(;
    name = "HT",
    decay_constant = TRITIUM_DECAY_CONSTANT,
    deposition_velocity = DepositionVelocity(0.04e-2, 0.05e-2),
)
