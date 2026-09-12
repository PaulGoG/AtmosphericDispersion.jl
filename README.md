# AtmosphericDispersion.jl

Gaussian-plume atmospheric dispersion of radioactive stack releases: dilution
factor, time-integrated air concentration, dry and wet ground deposition, and
resuspension, over Pasquill–Turner stability classes and a wind rose.

This branch is the rewrite and is not yet started. The solver core — plume rise,
Briggs dispersion parameters, building wake, deposition and resuspension — is
nuclide-agnostic; tritium (HTO/HT) is a parameterisation of it.

## Status

Nothing here yet beyond the licence. The work it will carry:

- Package scaffolding via `PkgTemplates.jl`, with the solver behind a typed,
  documented public interface and thin script entry points
- TOML-driven configuration, validated on load, replacing the hardcoded
  constants
- **Correcting the wind-rose sector convention.** The 2021 code numbers sectors
  counterclockwise from East, starting at a sector edge. The meteorological
  convention — which ADMS and the frequency tables of the governing norm both
  follow — centres sector 1 on North, increases clockwise, and denotes the
  direction the wind blows *from*. That is an axis rotation, a handedness flip
  and a half-sector binning offset, and it propagates into every long-duration
  result
- Removing the per-grid-point recomputation of building-equivalent geometry and
  the DataFrame mask lookups from the innermost loops
- Physics validation: dimensional analysis, limiting cases, and comparison
  against an independent implementation

## History

`original` holds the code exactly as submitted for the BSc thesis at the Faculty
of Physics, University of Bucharest, in June 2021, together with the thesis
itself. This branch shares no history with it.

## Licence

MIT, see `LICENSE`.
