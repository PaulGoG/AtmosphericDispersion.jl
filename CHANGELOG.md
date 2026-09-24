# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Four dispersion schemes**, selected by `[model] dispersion` and carried by
  `Site`. `hosker`, the normative's pairing of Briggs' open-country σ_y with
  Hosker's roughness-corrected σ_z, remains the default. `briggs_open_country`
  and `briggs_urban` take both parameters from Briggs (1973) as the Handbook on
  Atmospheric Diffusion prints them in Table 4.5, the urban set being his fit
  to the St. Louis experiment. `eimutis_konicek` is the fit of Eimutis and
  Konicek (1972) to the Pasquill–Gifford curves that NRC XOQDOQ and PAVAN
  evaluate, its 42 coefficients asserted against both codes' listings. The
  Handbook's `0.00015` inside the E–F urban σ_z is a misprint for Briggs'
  `0.0015`, checked in his own report and in EPA ISC3. `dispersion_parameters`
  returns both σ of a scheme; the existing `lateral_dispersion` and
  `vertical_dispersion` calls are unchanged.
- **Validity ranges.** Every scheme states the range of distance it was fitted
  over — 100 m to 10 km for the three Briggs-based schemes, 100 m to 100 km for
  Eimutis–Konicek — as `validity_range`, with `within_validity` for one
  distance. A configuration whose receptor grid reaches outside the band is
  loaded with a warning, refused, or accepted silently as
  `[model] extrapolation = "warn" | "refuse" | "allow"` directs, and
  `scripts/run.jl` reports how many radii extrapolate. The reference
  configuration's 20 km grid is twice the default band and says `allow`.
- `layer_mean_wind_speed`, the mean of the power-law profile over a layer in
  closed form, and `stable_rise_wind_speed`, the wind the stable final rise of
  a site is evaluated with.
- The validation figure `dispersion_schemes.png`: the four schemes side by side
  for class D, across their bands.
- **The building zones of IAEA SRS-19** §3.3 — displacement above `2.5 H_B`,
  wake beyond `2.5 √A_B`, cavity between — as `building_zone`, and the two
  cavity forms of §3.6: `dilution_cavity_wall`, `B₀/(u x²)` with `B₀ = 30` for
  a receptor on the wall the vent is in, and `dilution_cavity`,
  `1/(π u H_B K)` for one that is not.
- **Two more end-to-end benchmarks**, in `test/test_regulatory.jl`. IAEA SRS-19
  Tables I and II: all 77 cells of the diffusion factor and 106 of the 110 of
  the wake-corrected one reproduce to the printed figure, the four others
  sitting on a rounding boundary, and the three worked examples of Annex IV
  come out. NRC XOQDOQ Test Case 2, the code's own continuous elevated
  release — a momentum jet over rising terrain, summed over five wind classes
  and five stability classes: all 22 printed χ/Q reproduce to within 0.05 %,
  the rounding of the four printed figures.

### Changed

- **The stable final rise takes the mean wind over the rise.** The Handbook
  defines the `u` of Briggs' `2.6 [F/(uS)]^(1/3)`, its Eq. 2.19, as "an average
  value between the heights h_s and h_s + Δh"; the 2021 code, and XOQDOQ, took
  the wind at the release height. `Site` now solves the rise and the mean
  together, `[model] stable_rise_wind = "mean_over_rise"` by default, and
  `"release_height"` restores the old convention; the functions that take a
  scalar wind gain a `stable_wind` keyword. The reference atmosphere is stably
  stratified in every class, so every effective height falls by 1 to 3 m —
  class D at 2 km from 108.0 to 106.0 m, class F from 103.5 to 100.8 m, still
  above its 100 m lid — and the ground-level maxima rise by up to 8 % and move
  inward, class D from 1.51 × 10⁻⁶ s m⁻³ at 2.20 km to 1.58 × 10⁻⁶ at 2.15 km.
  The long-term χ/Q in the most exposed sector rises by 5 % at 1 km and by
  0.3 % at 20 km, where it stays 5.2 × 10⁻⁹ s m⁻³. The figures and the numbers
  in the documentation are regenerated.

## [0.2.0] - 2026-09-23

These change results, and several of them change the interface. Of them only
the mixing layer moves the reference case, by up to +23 % at 20 km: its stack is
buoyancy-dominated, so the plume-rise constants that changed are never reached;
it has no buildings; and its roughness class is not one of the two corrected.

### Added

- **A mixing layer.** The vertical profile between the ground and a capping
  inversion is HPA-RPD-058 §3.2.2.1 Eq. (3.4), summed to double precision —
  the image series while `Σ_z < 0.7 A`, its cosine series above, whose leading
  term is the uniform `1/A` of Eq. (3.5). Activity is conserved over
  `0 ≤ z ≤ A` to within 4 × 10⁻¹⁶ at every σ_z. `MixingLayer` holds one depth
  per stability class; the default is Table 3.5(a): A 1300, B 900, C 850,
  D 800, E 400, F 100 m. In a configuration it is the `[mixing_layer]` table:
  `scheme = "tabulated" | "uniform" | "custom" | "unbounded"`.
  - It raises ground-level χ/Q, never lowers it: 1.45× in class A at 20 km and
    2.23× at 50 km, 1.67× in class B at 50 km, and not at all in class D. The
    unstable classes reach the lid first because σ_z grows fastest there.
  - A plume above the lid is treated by a stated rule, `above_lid`.
    `rise_inhibited`, the default, holds it at the lid after NRPB-R157 §B2.3;
    it is continuous in the release height and is the case HPA-RPD-058
    Table 3.7 tabulates. `full_penetration` is EPA ISC3 vol. II §1.1.6.1: the
    plume leaves the layer and the ground sees nothing.
  - In the reference case class F rises to 103.5 m under a 100 m lid and is held
    there, which about doubles its ground-level value. The long-term χ/Q at
    20 km in the most exposed sector is 5.2 × 10⁻⁹ s m⁻³ against 4.2 × 10⁻⁹
    with no lid (undepleted).
- **`PrescribedPlume`**, a `Site` with its effective height, transport wind
  speed or dispersion parameters stated rather than derived. Published
  benchmarks give these as inputs, and a σ read off a graph is an input to the
  problem. The depletion integral refuses a prescribed σ_z, which holds at one
  distance only.
- **An end-to-end benchmark against HPA-RPD-058 Table 3.7** — 3 release heights
  × 6 stability categories × 8 distances of published plume-depletion fractions.
  143 of the 144 cells reproduce to within 0.028 against a printed precision of
  0.01, and every row is monotone. The exception is an error in the source:
  category F at 30 m reads 0.13 at 50 km and 0.19 at 100 km, a depletion factor
  that rises with distance; this package gives 0.0188, one decimal point from
  the printed value.
- **A second end-to-end benchmark, Turner Table 7-4** — the concentration
  profile with height from the ground to 450 m at 1 km, the only published table
  found that exercises both reflection terms at arbitrary receptor height. All
  16 rows reproduce to within 2.3 %, the residual being Turner's own
  three-significant-figure rounding.
- **Both species rows of the washout table.** NSR-23 Table 7 tabulates tritium
  and iodine separately from all other radionuclides, and only the first was
  carried. The row is a property of the species: `Nuclide.washout_species`,
  `[nuclide] washout_species`.
- **Snow washout for HTO**, `[model] washout = "hto"`, from Ogram (Ontario Hydro
  85-233-K, 1985) Eq. (38). Snow scavenging of tritiated water is isotopic
  exchange, not impaction. Ogram's values agree with NSR-23's own
  other-nuclides snow lower limit to within 16 % at all four rates. Rain is
  unaffected.
- **`ResuspensionModel`**, the two-exponential factor with a floor, carrying
  IAEA Safety Series 57 (the default) and Maxwell and Anspaugh (2011), selected
  by `[model] resuspension`; any other parameter set can be constructed.
- **Sources.** The tritium deposition velocities are those of the notes to
  NSR-23 Table 6, attributed there to Murphy, *Health Physics* **65**(6), 1993.
  The vertical dispersion scheme is Hosker's fit (IAEA-SM-181/19, 1974), printed
  in HPA-RPD-058 Table 3.3 and NRPB-R91 Table 3, and both tables are asserted
  entry by entry. The long-term sector constant is asserted against NRC
  Regulatory Guide 1.111 (2.032, sixteen sectors) and IAEA SRS-19 (1.5238,
  twelve). Safety Series 57's resuspension constants are asserted against the
  primary document.
- **A bibliography.** The sources the documentation and the docstrings cite
  are entries of `docs/src/refs.bib`, cited by key through DocumenterCitations
  and listed on a References page, with a DOI wherever one exists: Briggs
  (1969), Hanna, Briggs and Hosker (1982), Sagendorf, Goll and Sandusky (1982),
  Murphy (1993), Slinn (1977), Maxwell and Anspaugh (2011) and NUREG/CR-7270;
  the regulatory documents and laboratory reports by report number.

### Changed

- **Washout is one argument.** `WashoutEvent(; duration, precipitation, rate,
  model)` replaces the separate keywords on the dilution, depletion and
  deposition functions, and the matching fields of `RunConfiguration`. It is
  optional on the dilution factors (none means dry air) and required by the
  wet-deposition functions; passing one without a nuclide is an error.
- **Plume-rise constants default to Briggs, and the choice is configurable.**
  `RiseCoefficients` carries the three constants published schemes disagree on,
  selected by `[model] plume_rise`: `briggs` (the default), `xoqdoq` (stable
  coefficient 2.4), or `nsr23` (the 2021 values, `NSR23_RISE`).
  - The combined-law buoyancy denominator moves from 0.5 to 2β² = 0.72, so the
    law reduces to the two-thirds law as the momentum flux vanishes. It
    previously overshot by 13.6 %.
  - The neutral momentum rise moves from `1.5 w₀D/u` to Briggs' `3 w₀D/u`,
    Briggs (1969) Eq. 5.2 and EPA ISC3 Eq. (1-16). The 1.5 is the momentum term
    of Holland's formula (1953), Eq. (4.1) of Turner's Workbook.
- **The building-wake coefficient defaults to 1.0** rather than 1.5, which is
  IAEA SRS-19 Eq. (6) and the German AVV Eqs. (4.31)/(4.32). The 1.5 is that of
  CNCAN NSR-23, the Romanian normative the thesis followed;
  `NSR23_WAKE_COEFFICIENT` restores it.
- **Unknown configuration keys are an error** naming the key, where a misspelt
  optional key used to take its default in silence.
- `plume_rise(x, source, atmosphere, u, rise)` takes the rise coefficients and
  shares one implementation with the `Site` method.
- Julia 1.10, the LTS, is the floor. No `Manifest.toml` is tracked; every
  environment resolves against `[compat]`.
- The scripts and `docs/make.jl` activate their own environment, so
  `julia scripts/run.jl config/reference.toml` runs as written, without
  `--project`. The figure scripts take an optional output directory.
- `scripts/run.jl` evaluates the receptor grid of `[grid]`, `spacing` to
  `extent` in every sector, and reports the maximum over the radii beside the
  values at the extent, then the air concentration resuspended from the most
  exposed receptor's deposition under `[model] resuspension`. Both keys were
  validated on load and read by nothing.
- The comments of `config/reference.toml` state what each key is, its unit and
  its bounds. What the choices mean is in the configuration page of the manual.
- The narrative of the README moved into the manual: conventions, validation,
  and what the rewrite changed relative to the 2021 code.
- The σ_y citation named an NRC accession number that is not an NRC document.
  Corrected to Hanna, Briggs and Hosker, *Handbook on Atmospheric Diffusion*,
  DOE/TIC-11223 (1982), Table 4.5, which attributes ATDL Contribution No. 79.

### Fixed

- **Two coefficients of the roughness correction.** `F(z₀,x)` carried 1.58 and
  2.08 for z₀ = 0.01 m and 0.04 m where Hosker publishes 1.56 and 2.02, making
  σ_z 1.9 % too large over grassland and water and 3.6 % too large over arable
  land. These are the values NSR-23 Table 2 itself prints; the 2021 code carried
  them faithfully. The remaining 22 coefficients, and all 24 of the shape
  function, were already exact in both.

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
