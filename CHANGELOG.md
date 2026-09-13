# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0]

First public release. The solver is complete and validated; the API is not yet
stable, which is what the leading zero means.

### Added

- Sector geometry (`SectorGrid`) indexed clockwise from north and centred on the
  cardinal directions, and `WindRose`, which requires the wind-direction
  convention to be stated rather than assumed.
- Pasquill–Gifford stability classes, surface and roughness categories, and the
  coefficient tables of the governing normative as compile-time constants.
- Power-law wind profile, lateral and vertical dispersion parameters, plume rise
  by momentum and buoyancy, and the building-wake corrections.
- `Site`, which precomputes everything independent of receptor position.
- Three dilution regimes — instantaneous, extended and long term.
- Nuclides, plume depletion by decay, dry deposition and washout, and ground
  deposition with resuspension.
- TOML configuration validated key by key, with failures named by dotted path.
- Validation against the analytic invariants of the Gaussian plume and against
  the published Briggs (1973) parameterisations.

### Fixed, relative to the 2021 thesis code this derives from

- **Wind-direction convention.** The thesis code consumed a blowing-from wind
  rose as blowing-toward, which rotates the long-term dose field by half a turn
  and moves the most-exposed sector to the least-exposed one.
- **Sector binning.** Sectors began at an edge, putting every cardinal direction
  on a boundary where floating-point round-off chose the bin; two of sixteen
  sectors were unreachable and two received double weight.
- **Sector frame.** Sectors were indexed counterclockwise from east against a
  compass that runs clockwise from north — a reflection as well as a rotation.
- **Depletion composition.** Wet and dry depletion factors were added rather
  than multiplied, giving twice the undepleted concentration in the limit of no
  depletion at all.
- **Wet deposition normalisation.** `√2·π` where `√(2π)` was meant, a factor of
  `√π`.
- **Stability parameter.** The dry adiabatic lapse rate was computed with the
  specific heat of the emitted vapour rather than of the ambient air.
- Unstable atmospheres raised a domain error, since the stratification-limited
  plume-rise branches take fractional powers of a negative stability parameter.
- A plume trapped in a building cavity was transported at zero wind speed, which
  every dilution factor divides by.
