```@meta
CurrentModule = AtmosphericDispersion
```

# Configuration

A run is described by a TOML file, never by editing source. The loader checks
presence, type, enumerated choice and numerical bound, rejects any key it does
not read, and attaches every failure to the key it failed on by its full dotted
path.

```julia-repl
julia> load_configuration("broken.toml")
ERROR: ConfigurationError at `grid.spacing`: must be smaller than grid.extent (20000.0 m), got 20000.0

julia> load_configuration("misspelt.toml")
ERROR: ConfigurationError at `model.mixing_layr`: unknown key; expected one of plume_rise, resuspension, mixing_layer, washout
```

An unread key is an error because it cannot be told from a misspelt one, and a
misspelt optional key would otherwise take its default without a word.

Validation runs at load and constructs the solver types immediately, so
anything those constructors reject in turn — frequencies that do not sum to
one, a negative stack height — fails at the same point rather than deep into a
run.

`config/reference.toml` is the annotated template: every key, its unit and its
bounds. The comments there are limited to that; what the choices mean is below.

## The `[model]` table

Where published schemes disagree, the choice is made here. Every key is
optional and defaults to the first value listed.

| Key | Value | Meaning |
|---|---|---|
| `plume_rise` | `briggs` | [`BRIGGS_RISE`](@ref): `2β² = 0.72`, `3 w₀D/u`, `2.6 [F/(uS)]^(1/3)` |
| | `xoqdoq` | [`XOQDOQ_RISE`](@ref): as `briggs`, with the stable coefficient 2.4 of NUREG/CR-2919 [Sagendorf1982](@cite) |
| | `nsr23` | [`NSR23_RISE`](@ref): 0.5 and 1.5, the constants of CNCAN NSR-23 [CNCAN2004](@cite) as the 2021 code carried them |
| `stable_rise_wind` | `mean_over_rise` | `WIND_MEAN_OVER_RISE`: the `u` of the stable final rise is the mean of the profile between the release height and the top of the rise, as the Handbook [Hanna1982](@cite) defines it for its Eq. 2.19 |
| | `release_height` | `WIND_AT_RELEASE_HEIGHT`: the transport wind at the release height, as XOQDOQ and the 2021 code took it |
| `dispersion` | `hosker` | `DISPERSION_HOSKER`: Briggs' open-country `σ_y` with Hosker's roughness-corrected `σ_z`, the normative's pairing |
| | `briggs_open_country` | `DISPERSION_BRIGGS_OPEN_COUNTRY`: both parameters from [Briggs1973](@citet), Handbook Table 4.5 |
| | `briggs_urban` | `DISPERSION_BRIGGS_URBAN`: Briggs' urban set, fitted to the St. Louis experiment of [McElroyPooler1968](@citet) |
| | `eimutis_konicek` | `DISPERSION_EIMUTIS_KONICEK`: the fit of [EimutisKonicek1972](@citet) to the Pasquill–Gifford curves, as NRC XOQDOQ and PAVAN [Bander1982](@cite) evaluate it |
| `extrapolation` | `warn` | a receptor grid reaching outside the scheme's [`validity_range`](@ref) is loaded with a warning naming the offending end |
| | `refuse` | the same grid is refused, on `grid.spacing` or `grid.extent` |
| | `allow` | the grid is accepted silently; the reference configuration says this, its 20 km extent being twice the 10 km band of the default scheme |
| `resuspension` | `iaea_ss57` | IAEA Safety Series 57 [IAEA1982](@cite), Eq. (3.14A) |
| | `maxwell_anspaugh` | [MaxwellAnspaugh2011](@citet); NUREG/CR-7270 [Bixler2022](@cite) |
| `washout` | `normative` | NSR-23 Table 7 throughout |
| | `hto` | snow scavenging from [Ogram1985](@citet) Eq. (38), for tritiated water; rain unchanged |

See [`RiseCoefficients`](@ref), [`StableRiseWind`](@ref),
[`DispersionScheme`](@ref), [`ResuspensionModel`](@ref) and
[`WashoutModel`](@ref) for the physics behind each. The validity range of a
scheme is the range of distance its parameters were fitted over; every
Briggs-based scheme is quoted for 100 m to 10 km, the Eimutis–Konicek fit for
100 m to 100 km. Nothing in the evaluation refuses a distance outside it; the
policy applies to the receptor grid at load, and `scripts/run.jl` reports how
many radii extrapolate.

## The `[mixing_layer]` table

| Key | Value | Meaning |
|---|---|---|
| `scheme` | `tabulated` | HPA-RPD-058 [SmithSimmonds2009](@cite) Table 3.5(a) by class: A 1300, B 900, C 850, D 800, E 400, F 100 m |
| | `uniform` | one depth for every class, given as `uniform_depth` in metres; HPA-RPD-058 §3.2.2.2.3 recommends 800 |
| | `custom` | six `depths` in metres, class A to F; `inf` leaves a class unbounded |
| | `unbounded` | no lid |
| `above_lid` | `rise_inhibited` | a plume above the lid is held at it, NRPB-R157 [Jones1983](@cite) §B2.3 |
| | `full_penetration` | a plume above the lid leaves the layer and gives no ground-level concentration, EPA ISC3 vol. II [EPA1995](@cite) §1.1.6.1 |

`uniform_depth` is required by, and only legal with, `scheme = "uniform"`, and
`depths` likewise with `scheme = "custom"`. See [`MixingLayer`](@ref) and
[`LidRule`](@ref).

## The `[nuclide]` table

`deposition_velocity_low` feeds the depletion of the airborne plume and
`deposition_velocity_high` the ground deposition, each conservative for the
endpoint it feeds; see [`DepositionVelocity`](@ref). The reference values,
0.4–0.8 × 10⁻² m/s for HTO, are those of the notes to NSR-23 Table 6,
attributed there to [Murphy1993](@citet). IAEA SRS-19 [IAEA2001](@cite) §3.9,
EUR 15760 [Simmonds1995](@cite) §3.2 and HPA-RPD-058 §3.2.2.3 instead assign tritium no
dry deposition at all; set both velocities to zero to follow them.

`washout_species` selects the row of NSR-23 Table 7: `tritium_iodine`, the
default, or `other` for every other radionuclide. See [`WashoutSpecies`](@ref).

## The `[buildings]` table

`wake_coefficient` is `C` in `√(σ² + C A/π)`. The default of 1 is IAEA SRS-19
Eq. (6) and AVV [AVV2012](@cite) Eqs. (4.31)/(4.32); NSR-23 uses 1.5
([`NSR23_WAKE_COEFFICIENT`](@ref)). Zero disables the correction.

## The `[wind_rose]` table

`convention` is required and has no default: `blowing_from` for a measured rose
or an ADMS `.met` file, `blowing_toward` for frequencies already expressed as
transport towards a sector. See [`WindRose`](@ref).

!!! warning "The reference frequencies are illustrative"
    The rose in `config/reference.toml` is not a measured rose for any one
    site. I assembled it in 2021 from meteorological data for the Pitești fuel
    plant and the Cernavodă NPP, and it is unusually flat for a real site.
    Replace it before using the configuration for an assessment.

## Reference

```@autodocs
Modules = [AtmosphericDispersion]
Pages = ["config.jl"]
```
