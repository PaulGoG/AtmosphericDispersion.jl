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
| | `xoqdoq` | [`XOQDOQ_RISE`](@ref): as `briggs`, with the stable coefficient 2.4 of NUREG/CR-2919 |
| | `nsr23` | [`NSR23_RISE`](@ref): 0.5 and 1.5, the constants of CNCAN NSR-23 as the 2021 code carried them |
| `resuspension` | `iaea_ss57` | IAEA Safety Series 57 (1982), Eq. (3.14A) |
| | `maxwell_anspaugh` | Maxwell and Anspaugh, *Health Physics* **101** (2011); NUREG/CR-7270 |
| `washout` | `normative` | NSR-23 Table 7 throughout |
| | `hto` | snow scavenging from Ogram (1985) Eq. (38), for tritiated water; rain unchanged |

See [`RiseCoefficients`](@ref), [`ResuspensionModel`](@ref) and
[`WashoutModel`](@ref) for the physics behind each.

## The `[mixing_layer]` table

| Key | Value | Meaning |
|---|---|---|
| `scheme` | `tabulated` | HPA-RPD-058 Table 3.5(a) by class: A 1300, B 900, C 850, D 800, E 400, F 100 m |
| | `uniform` | one depth for every class, given as `uniform_depth` in metres; HPA-RPD-058 §3.2.2.2.3 recommends 800 |
| | `custom` | six `depths` in metres, class A to F; `inf` leaves a class unbounded |
| | `unbounded` | no lid |
| `above_lid` | `rise_inhibited` | a plume above the lid is held at it, NRPB-R157 §B2.3 |
| | `full_penetration` | a plume above the lid leaves the layer and gives no ground-level concentration, EPA ISC3 vol. II §1.1.6.1 |

`uniform_depth` is required by, and only legal with, `scheme = "uniform"`, and
`depths` likewise with `scheme = "custom"`. See [`MixingLayer`](@ref) and
[`LidRule`](@ref).

## The `[nuclide]` table

`deposition_velocity_low` feeds the depletion of the airborne plume and
`deposition_velocity_high` the ground deposition, each conservative for the
endpoint it feeds; see [`DepositionVelocity`](@ref). The reference values,
0.4–0.8 × 10⁻² m/s for HTO, are those of the notes to NSR-23 Table 6. IAEA
SRS-19 §3.9, EUR 15760 §3.2 and HPA-RPD-058 §3.2.2.3 instead assign tritium no
dry deposition at all; set both velocities to zero to follow them.

`washout_species` selects the row of NSR-23 Table 7: `tritium_iodine`, the
default, or `other` for every other radionuclide. See [`WashoutSpecies`](@ref).

## The `[buildings]` table

`wake_coefficient` is `C` in `√(σ² + C A/π)`. The default of 1 is IAEA SRS-19
Eq. (6) and AVV Eqs. (4.31)/(4.32); NSR-23 uses 1.5
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
