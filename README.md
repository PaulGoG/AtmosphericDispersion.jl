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
├── src/
│   ├── AtmosphericDispersion.jl   module, includes and exports
│   ├── sectors.jl          wind-rose sector geometry
│   ├── stability.jl        Pasquill–Gifford stability classes
│   └── windrose.jl         directional and stability joint frequencies
├── test/
│   ├── Project.toml
│   ├── activate.jl
│   └── runtests.jl         unit tests + Aqua, JET, ExplicitImports
├── docs/
│   ├── Project.toml
│   ├── activate.jl
│   └── make.jl             Documenter build
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
julia --project -e 'using Pkg; Pkg.test()'          # test suite and static QA
julia --project=docs docs/make.jl                    # build the documentation
julia -e 'using JuliaFormatter; format(".")'         # apply the committed style
```

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

## Status

Foundations in place and under test; the solver itself is not yet written.

Done:

- Sector geometry (`SectorGrid`), Pasquill classes, and the wind rose with its
  explicit direction convention
- Static QA in the suite: Aqua, JET, ExplicitImports

Next:

- TOML-driven configuration, validated on load, replacing hardcoded constants
- Plume rise, Briggs dispersion parameters and building-wake corrections
- Depletion, dry and wet deposition, resuspension
- The three dilution regimes behind a typed public interface, with thin script
  entry points
- Removal of the per-grid-point recomputation of building-equivalent geometry
  and of the DataFrame mask lookups in the innermost loops
- Physics validation: dimensional analysis, limiting cases, and comparison
  against an independent implementation

## History

`original` holds the code exactly as submitted for the BSc thesis at the Faculty
of Physics, University of Bucharest, in June 2021, together with the thesis
itself. This branch shares no history with it.

## Licence

MIT, see `LICENSE`.
