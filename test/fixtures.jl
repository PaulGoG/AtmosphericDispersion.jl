# The reference case the suite is evaluated on: the stack and atmosphere of the
# 2021 run. The constructor functions take the same keywords as the package
# types, so a test states only what it changes.

reference_stack(;
    height = 50.3,
    diameter = 2.33,
    exit_velocity = 10.0,
    exit_density = 0.6,
    exit_temperature = 324.0,
) = StackSource(; height, diameter, exit_velocity, exit_density, exit_temperature)

reference_atmosphere(;
    reference_speed = 4.0,
    temperature = 287.0,
    density = 1.2,
    lapse_rate = 2e-2,
    surface = SURFACE_AGRICULTURAL,
    roughness = ROUGHNESS_PASTURE,
) = Atmosphere(; reference_speed, temperature, density, lapse_rate, surface, roughness)

const REFERENCE_STACK = reference_stack()
const REFERENCE_AIR = reference_atmosphere()
const REFERENCE_SITE = Site(; source = REFERENCE_STACK, atmosphere = REFERENCE_AIR)

# The same site with no mixing lid, for the statements that hold only for the
# unbounded Gaussian.
const UNBOUNDED_SITE =
    Site(; source = REFERENCE_STACK, atmosphere = REFERENCE_AIR, mixing = MIXING_UNBOUNDED)

# A species that neither decays nor deposits, so every depletion factor is one.
const INERT = Nuclide(;
    name = "inert",
    decay_constant = 0.0,
    deposition_velocity = DepositionVelocity(0.0, 0.0),
)

# Stability fractions A to F of the 2021 input table, as given there: they sum to
# one only to the five figures they were recorded at.
const STABILITY_2021 = (0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158)

# The same, normalised.
normalised_stability() = collect(STABILITY_2021) ./ sum(STABILITY_2021)
