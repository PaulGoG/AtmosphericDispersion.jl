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
│   └── AtmosphericDispersion.jl
├── test/
│   ├── Project.toml
│   ├── activate.jl
│   └── runtests.jl         unit tests + Aqua static QA
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

## Status

Scaffolding only. The solver is not yet written. What this branch will carry:

- TOML-driven configuration, validated on load, replacing hardcoded constants
- The solver behind a typed, documented public interface, with thin script
  entry points
- **A corrected wind-rose sector convention.** The 2021 code numbers sectors
  counterclockwise from East, starting at a sector edge. The meteorological
  convention — which ADMS and the frequency tables of the governing norm both
  follow — centres sector 1 on North, increases clockwise, and denotes the
  direction the wind blows *from*. That is an axis rotation, a handedness flip
  and a half-sector binning offset, and it propagates into every long-duration
  result
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
