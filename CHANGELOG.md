# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **An end-to-end benchmark against HPA-RPD-058 Table 3.7** — 3 release heights
  × 7 stability categories × 8 distances of published plume-depletion fractions.
  143 of the 144 cells reproduce to within 0.028 against a printed precision of
  0.01, and every row is monotone. The exception is an error in the source:
  category F at 30 m reads 0.13 at 50 km and 0.19 at 100 km, a depletion factor
  that rises with distance; this package gives 0.0188, one decimal point from
  the printed value. The test asserts both the defect and the computed value.
- **A second end-to-end benchmark, Turner Table 7-4** — the concentration
  profile with height from the ground to 450 m at 1 km, the only published table
  found that exercises both reflection terms at arbitrary receptor height, where
  everything else evaluates at z = 0 and the two coincide. All 16 rows reproduce
  to within 2.4 %, the residual being Turner's own three-significant-figure
  rounding. His printed prefactor, 3.5 × 10⁻⁵ g/m³, is a typo: 151/(2π·157·110·4)
  is 3.479 × 10⁻⁴, which is what his own table uses.
- **`Site` takes `fixed_lateral` and `fixed_vertical`** as well, overriding the
  two dispersion parameters. Turner's problems print σ read off his figures, and
  a σ read off a graph is an input to the problem.
- **`Site` takes `fixed_height` and `fixed_wind`**, overriding the effective
  release height and transport wind speed. Published benchmarks state both as
  inputs rather than deriving them, so reproducing one requires setting them.
- **A mixing layer.** `MixingLayer`, selected by `[model] mixing_layer`,
  implements HPA-RPD-058 §3.2.2.1 Eq. (3.4) — the image sum over virtual sources
  at `2sA ± h_e`, truncated at `|s| = 1` — going over to Eq. (3.5), a uniform
  `1/A` profile, once σ_z reaches the depth. Depths are Table 3.5(a): A 1300,
  B 900, C 850, D 800, E 400, F 100 m, attributed there to Clarke (1979) and
  Jones (1980). `unbounded` restores the previous behaviour and `uniform_800`
  takes the single depth the report recommends for all conditions.
  - It raises ground-level χ/Q, never lowers it: 1.44× in class A at 20 km and
    2.23× at 50 km, 1.67× in class B at 50 km, and essentially not at all in
    class D. The unstable classes reach the lid first because σ_z grows fastest
    there, despite their deeper lids.
  - Activity is confined to `0 ≤ z ≤ A` and conserved exactly over it.
  - A release strictly above the lid is not trapped by it: class F's tabulated
    depth is 100 m against a 103 m effective release, so the plume starts above
    the inversion and is decoupled until it breaks — fumigation, which this
    package does not model. A release exactly at the lid is still capped, which
    Table 3.7 decided rather than taste.

### Changed

- **The building-wake coefficient defaults to 1.0** rather than 1.5, which is
  IAEA SRS-19 Eq. (6) and the German AVV Eqs. (4.31)/(4.32). The 1.5 came from
  CNCAN NSR-23, the Romanian normative the thesis followed, and no source
  outside it was found; `NORMATIVE_WAKE_COEFFICIENT` restores it.
- **Both species rows of the washout table.** NSR-23 Table 7 tabulates tritium
  and iodine separately from all other radionuclides, and only the first was
  carried; `WashoutSpecies` selects between them. The other-nuclides snow values
  run three to four orders of magnitude *above* their own rain values, opposite
  to the tritium row's suppression.
- **Snow washout for HTO is selectable**, `[model] washout = "hto"`, using Ogram
  (Ontario Hydro 85-233-K, 1985) Eq. (38). Snow scavenging of tritiated water is
  isotopic exchange, not impaction. Ogram's measurement agrees with NSR-23's own
  other-nuclides snow lower limit to within 16 % at all four rates, which is the
  evidence that the tritium row's snow column is the anomaly. Rain is unaffected.
- **Plume-rise constants now default to Briggs, and the choice is
  configurable.** `RiseCoefficients` carries the three constants published
  schemes disagree on, selected by `[model] plume_rise`: `briggs` (the default),
  `xoqdoq` (stable coefficient 2.4), or `thesis_2021` (the 2021 values).
  - The combined-law buoyancy denominator moves from 0.5 to **2β² = 0.72**, so
    the law reduces to the two-thirds law as the momentum flux vanishes. It
    previously overshot by 13.6 %.
  - The neutral momentum rise moves from `1.5 w₀D/u` to Briggs' **`3 w₀D/u`**,
    per Briggs (1969) Eq. 5.2 and EPA ISC3 Eq. (1-16). No source was found for
    1.5.
  - These change results. `plume_rise = "thesis_2021"` reproduces the old ones,
    and the fidelity tests use it.
- `resuspension_factor` takes a `ResuspensionModel`, selected by
  `[model] resuspension`: `iaea_ss57` (the default) or `maxwell_anspaugh`.
- `Site` gains a `rise` field; `RunConfiguration` gains `resuspension`.

### Fixed

- **Two coefficients of the roughness correction.** `F(z₀,x)` carried 1.58 and
  2.08 for z₀ = 0.01 m and 0.04 m where Hosker publishes 1.56 and 2.02, making
  σ_z 1.9 % too large over grassland and water and 3.6 % too large over arable
  land. These are the values **NSR-23 Table 2 itself prints** — the 2021 code
  carried them faithfully, and the normative attributes its table to
  CAN/CSA-N288.2-M91 and UNSCEAR 2000, so the corruption entered in that chain.
  The remaining 22 coefficients, and all 24 of the shape function, were already
  exact in both.

### Added

- Safety Series 57's resuspension constants and Table II washout constants, and
  NSR-23's own tables, are now asserted against the primary documents rather
  than against secondary reporting of them.
- The vertical dispersion scheme is identified: it is Hosker's fit
  (IAEA-SM-181/19, 1974) to F.B. Smith (1972) and Briggs (1973), printed in
  HPA-RPD-058 Table 3.3 and NRPB-R91 Table 3. Both tables are now asserted entry
  by entry against those printings rather than sampled.
- The long-term sector constant is asserted against its published value in NRC
  Regulatory Guide 1.111 (2.032, sixteen sectors) and IAEA SRS-19 (1.5238,
  twelve), for every stability class and distance.
- A test recording the one known deviation from Briggs: the combined
  momentum-and-buoyancy rise overshoots the two-thirds law by 13.6 % as the
  momentum flux vanishes, because its buoyancy denominator is 0.5 where Briggs
  gives 2β² = 0.72. The constant is left as the thesis set it; the test pins the
  size of the discrepancy so it cannot drift unnoticed.

### Changed

- The σ_y citation named an NRC accession number that is not an NRC document.
  Corrected to Hanna, Briggs and Hosker, *Handbook on Atmospheric Diffusion*,
  DOE/TIC-11223 (1982), Table 4.5, which attributes ATDL Contribution No. 79.

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
