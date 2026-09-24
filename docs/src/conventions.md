# Conventions

## The wind-direction convention

A wind direction carries a convention that its magnitude does not reveal, and
getting it wrong rotates a long-term dispersion field by half a turn while
leaving every magnitude, sum rule and unit intact. Three things had to be fixed
relative to the 2021 code, and they are independent of one another.

**Sense.** ADMS, and every measured wind rose, gives the direction the wind
blows *from*. The long-term sector-averaged formula weights by the frequency of
wind blowing *towards* the receptor's sector. `WindRose` therefore takes the
convention as a required argument — `BlowingFrom()` or `BlowingToward()` — and
stores blowing-towards internally, so the formulae cannot be fed the wrong
sense. There is no default.

**Frame.** The 2021 code indexed sectors from `atan(y/x)`, counterclockwise from
east. The compass runs clockwise from north. Mapping a cardinal table onto those
indices is a reflection as well as a rotation, so a table indexed N, NNE, NE, …
could not be fed in as k = 1, 2, 3, … under any sign convention.

**Binning.** The 2021 sectors began at a sector edge, which put every cardinal
direction exactly on a boundary, where `floor` on a ratio that came out one ulp
low decided the bin. ESE and SE collapsed into one sector, SSW and SW into
another, and two of the sixteen sectors became unreachable. Sectors here are
centred on the cardinal directions and binned by rounding to the nearest centre,
so a cardinal direction sits as far from a boundary as it can and bins exactly.

The thesis anticipated this last hazard and `Circle_sectors.jl` was written to
check it, but that script sampled a *random* point on the circle, which almost
surely never lands on a boundary — so the test could not detect the failure it
was written for. All sixteen cardinal directions, and both sides of every
boundary, are now asserted explicitly.

## What the convention costs

Reading the 2021 frequency table under the two conventions, everything else
held fixed, the long-term dilution factor at 1 km differs by up to a factor of
**1.63**, and the most-exposed sector moves from **N** to **S** — exactly
opposite, as it must. For a dose assessment the sector that matters is the
most-exposed one, so this is not a refinement: it points the assessment at the
wrong side of the site.

**I misread the convention in 2021.** The table is therefore a *blowing-from*
rose, which is what every published wind
rose is, and the 2021 code consumed it as *blowing-toward*. The 2021 dose field
is rotated by half a turn, and its most-exposed sector is the least-exposed one.
`config/reference.toml` declares `blowing_from` accordingly.

The provenance of the numbers is still open: as far as I recall they are either
placeholder values or real meteorological data for the Pitești fuel plant or the
Cernavodă NPP site. Read in standard cardinal order the rose is
north-dominated — 30.6 % of the time in the N quadrant against 19.8 % in the S,
strongest from N, NNE and E, weakest from SSW and SW — which is the right shape
for Dobrogea, where the crivăț blows from the north-east. That is consistent
with the Cernavodă possibility but does not establish it: the
rose is also unusually flat for a real site, only 1.7:1 between its strongest
and weakest sectors, where measured roses are typically more peaked. Treat the
absolute frequencies as unverified; the *convention* is not.

## Depletion composes by multiplication

A second defect, independent of the wind rose and provable the same way — by a
limit. The 2021 code combined the wet and dry depletion factors additively,

    χ = χ/Q · Q · DEC · (DEP_w + DEP_d)

Switch the rain off and set the deposition velocity to zero, so nothing is
removed at all: each factor tends to one, their sum tends to two, and the
concentration comes out twice the undepleted value. Surviving fractions
multiply — the processes act in sequence on the same material.

The long-term dry factor had the same shape of error at larger scale. It summed
six exponentials, one per stability class, with the frequencies inside the
exponent and no weighting outside, so in the no-deposition limit it tended to
six rather than one. Depletion is class-dependent through both the transport
speed and the vertical dispersion, so it belongs inside the class sum, and that
is where it now sits.

Both limits are asserted in the suite.

A third defect of the same kind sits in the wet deposition. Integrating the
Gaussian plume over the whole vertical column leaves a normalisation of
`√(2π) Σ_y u`; the 2021 code wrote `√2 π Σ_y u`. The ratio of the two is
exactly `√π`, the signature of `√(2π)` mistyped as `√2·π`, and it understated
wet deposition by a factor of 1.772.
