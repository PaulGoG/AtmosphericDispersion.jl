# AtmosphericDispersion.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PaulGoG.github.io/AtmosphericDispersion.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PaulGoG.github.io/AtmosphericDispersion.jl/dev/)
[![Build Status](https://github.com/PaulGoG/AtmosphericDispersion.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PaulGoG/AtmosphericDispersion.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/PaulGoG/AtmosphericDispersion.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PaulGoG/AtmosphericDispersion.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Gaussian-plume atmospheric dispersion of radioactive stack releases: dilution
factor, time-integrated air concentration, dry and wet ground deposition, and
resuspension, over Pasquill–Turner stability classes and a wind rose.

The solver core — plume rise, Briggs dispersion parameters, building wake,
deposition and resuspension — is nuclide-agnostic. Tritium (HTO/HT) release from
a CANDU stack is a parameterisation of it, and is the case the original work
addressed.

## Layout

```
.
├── Project.toml            package manifest and [compat]
├── activate.jl             activates and instantiates the root environment
├── config/
│   └── reference.toml      the CANDU tritium case, with bounds in the comments
├── scripts/
│   └── run.jl              long-term field for a configured run
├── src/
│   ├── AtmosphericDispersion.jl   module, includes and exports
│   ├── sectors.jl          wind-rose sector geometry
│   ├── stability.jl        Pasquill–Gifford stability classes
│   ├── windrose.jl         directional and stability joint frequencies
│   ├── surfaces.jl         surface, roughness and precipitation categories
│   ├── tables.jl           tabulated coefficients of the normative
│   ├── dispersion.jl       wind profile and dispersion parameters
│   ├── source.jl           stack, ambient state, fluxes, stability parameter
│   ├── plumerise.jl        momentum and buoyancy rise
│   ├── buildings.jl        building envelope and wake broadening
│   ├── site.jl             release height and the per-run precomputation
│   ├── dilution.jl         the three dilution regimes
│   ├── nuclides.jl         species: decay constant and deposition velocity
│   ├── depletion.jl        decay, dry deposition and washout
│   ├── deposition.jl       ground deposition and resuspension
│   └── config.jl           TOML loader, validated key by key
├── test/
│   ├── Project.toml
│   ├── activate.jl
│   └── runtests.jl         unit tests + Aqua, JET, ExplicitImports
├── docs/
│   ├── Project.toml
│   ├── activate.jl
│   ├── make.jl             Documenter build
│   └── src/                overview and four API pages
└── .github/workflows/      CI, CompatHelper, TagBot
```

## Environments

Each environment carries its own activation script, which activates and
instantiates it silently. Expect the first run of each to be slow.

```bash
julia activate.jl            # root
julia test/activate.jl       # test environment
julia docs/activate.jl       # documentation environment
```

## Entry points

```bash
julia --project scripts/run.jl config/reference.toml  # long-term field for a run
julia --project -e 'using Pkg; Pkg.test()'            # test suite and static QA
julia --project=docs docs/make.jl                      # build the documentation
julia -e 'using JuliaFormatter; format(".")'           # apply the committed style
```

A run is described by a TOML file, never by editing source. The loader checks
presence, type, enumerated choice and numerical bound, and names the offending
key by its full dotted path — `ConfigurationError at \`grid.spacing\`: must be
smaller than grid.extent (20000.0 m), got 20000.0` — so a rejected file says
what to change.

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

The thesis anticipated this last hazard and `SectoareCerc.jl` was written to
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

**The author has since confirmed he misunderstood the convention at the time.**
The table is therefore a *blowing-from* rose, which is what every published wind
rose is, and the 2021 code consumed it as *blowing-toward*. The 2021 dose field
is rotated by half a turn, and its most-exposed sector is the least-exposed one.
`config/reference.toml` declares `blowing_from` accordingly.

The provenance of the numbers is still open: his recollection is that they are
either placeholder values or real meteorological data for the Pitești fuel plant
or the Cernavodă NPP site. Read in standard cardinal order the rose is
north-dominated — 30.6 % of the time in the N quadrant against 19.8 % in the S,
strongest from N, NNE and E, weakest from SSW and SW — which is the right shape
for Dobrogea, where the crivăț blows from the north-east. That is consistent
with the Cernavodă branch of his recollection but does not establish it: the
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

## Status

The solver is complete and configuration-driven. Documentation and validation
against an independent implementation remain.

Done:

- Sector geometry, Pasquill classes, and the wind rose with its explicit
  direction convention
- Surface and roughness categories; the coefficient tables as compile-time
  constants rather than DataFrame lookups in the inner loop
- Wind profile and dispersion parameters; plume rise; building wake
- `Site`, which precomputes everything independent of receptor position
- The three dilution regimes: instantaneous, extended, long term
- Nuclides, and depletion by decay, dry deposition and washout
- Dry and wet ground deposition, and resuspension
- TOML configuration validated key by key, and a script entry point over it
- Static QA in the suite: Aqua, JET, ExplicitImports. 2920 tests, including
  exact agreement with the 2021 formulae wherever they were evaluable
- A Documenter site that builds clean, doctests included

Next:

- Comparison against an independent implementation
- Establish the provenance of the 2021 wind rose, and settle its convention

## History

`original` holds the code exactly as submitted for the BSc thesis at the Faculty
of Physics, University of Bucharest, in June 2021, together with the thesis
itself. This branch shares no history with it.

## Licence

MIT, see `LICENSE`.
