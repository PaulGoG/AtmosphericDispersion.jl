# Validation

The solver is checked against the analytic properties of the Gaussian plume, not
only against itself. These are the statements that fail first if a
normalisation, a reflection term or a sector width is wrong.

| Check | Result |
|---|---|
| Crosswind integral equals `2exp(−H²/2Σ_z²)/(√(2π) Σ_z u)` | exact to 1 part in 10⁸ |
| Mass conservation, `u ∫∫ χ/Q dy dz = 1` | **1.000000000** |
| Sector average equals the crosswind integral over the arc | exact to 1 part in 10¹⁰ |
| Ground-level maximum at `Σ_z = H/√2` | exact under its own assumptions |

Mass conservation is the strongest of them: it says the released activity
crossing any downwind plane, carried at the transport speed, is the whole
release — which exercises the 2πΣ_yΣ_z normalisation, the ground reflection and
the transport speed together.

## Against the published literature

The governing normative is **CNCAN NSR-23** [CNCAN2004](@cite), the Romanian nuclear regulator's
norm the 2021 thesis was written against. Every parameterisation it hands down
turns out to be a standard published scheme, and where the two disagree this
package now follows the international source and keeps the normative's value
available. Where they are, the identity is asserted rather than
described.

| | Published form | Agreement |
|---|---|---|
| σ_y | [Briggs1973](@citet) open country, `a x(1+10⁻⁴x)^(−1/2)`, a = 0.22…0.04 | **exact**, all six classes |
| σ_z shape `g(x)` | [Hosker1974](@citet), `a x^b/(1 + c x^d)` | **exact**, all 24 coefficients |
| σ_z roughness `F(z₀,x)` | Hosker (1974), both branches | **exact**, all 24 coefficients — after correcting two |
| Long-term sector constant | NRC RG 1.111 [NRC1977](@cite) `2.032`, IAEA SRS-19 [IAEA2001](@cite) `1.5238` | **exact to the published rounding** |
| Distance to final rise | Briggs, as in [Hanna1982](@citet), `x_f = 14F^(5/8)`, `34F^(2/5)` | **exact** |
| Neutral final buoyant rise | Briggs `21.4F^(3/4)/u`, `38.7F^(3/5)/u` | **exact up to the literature's own rounding** |
| Stable final rise | Briggs `2.6[F/(us)]^(1/3)` | **exact** |
| Transitional rise | the two-thirds law `1.6F^(1/3)x^(2/3)/u` | **exact** |
| Combined rise | Briggs, `3F_m x/β_j²u² + 3F x²/2β²u³` | momentum half exact, **buoyancy half 13.6 % high** |

### The vertical dispersion scheme, and two wrong digits

`σ_z = g(x)F(z₀,x)` is Hosker's analytic fit [Hosker1974](@cite) to the
schemes of F.B. Smith [Smith1973](@cite) and Briggs [Briggs1973](@cite),
published as IAEA-SM-181/19 in *Physical Behaviour of Radioactive Contaminants
in the Atmosphere*. The normative reproduces it without attribution. It is
printed in full in **HPA-RPD-058** [SmithSimmonds2009](@cite) Table 3.3 and
**NRPB-R91** [Clarke1979](@cite) Table 3, and both printings agree digit for
digit.

Asserting the whole table against those printings found **two wrong
coefficients** in the roughness correction:

| z₀ | NSR-23 Table 2 | Hosker | effect on σ_z |
|---|---|---|---|
| 0.01 m, grassland and water | 1.58 | **1.56** | 1.9 % too large |
| 0.04 m, arable | 2.08 | **2.02** | 3.6 % too large |

**The error is the normative's, not the thesis's.** NSR-23 prints 1.58 and 2.08,
and the 2021 code carried them faithfully. NSR-23 attributes its Table 2 to
CAN/CSA-N288.2-M91 and UNSCEAR 2000, so the corruption entered somewhere in that
chain rather than at either end. This package follows Hosker.

The other 22 coefficients, and all 24 of the shape function, match exactly in
both NSR-23 and Hosker. This is the argument for asserting tables of constants
entry by entry rather than sampling them: nothing else in the suite could have
caught it, because the error is in the data, not the algebra.

### The long-term equation is a regulatory one

`dilution_long_term` is not a bespoke form. The same equation, constant
included, is stated by **NRC Regulatory Guide 1.111** Rev. 1 [NRC1977](@cite)
Eq. (3) and its implementation **XOQDOQ** [Sagendorf1982](@cite) Eq. (1), by
**IAEA Safety Reports Series No. 19** [IAEA2001](@cite) Eq. (V-2), and by the
German **AVV zu §47 StrlSchV** [AVV2012](@cite) Eq. (4.4). RG 1.111 writes the sector constant as `2.032` and
says in words that it is `√(2/π)` divided by a 22.5° sector in radians; SRS-19
works in twelve sectors, where the same constant is `1.5238`. Both are asserted,
for every class and distance.

`dilution_instantaneous` is likewise SRS-19 Eq. (V-1) literally.

### An end-to-end benchmark

Everything above checks one piece at a time. **HPA-RPD-058 Table 3.7** checks
them together: *"Fractions of material remaining in the plume due to dry
deposition for a deposition velocity of 10⁻² m s⁻¹"*, three effective release
heights × seven stability categories × eight distances out to 100 km, with the
stack-height wind speed printed for each row.

Reproducing it exercises σ_z, the mixing layer and the depletion integral at
once. **143 of the 144 cells agree to within 0.028**, the table's own printed
precision being 0.01, and every row is monotone.

The one exception is not the package. Category F at 30 m reads **0.13 at 50 km
and 0.19 at 100 km** — a depletion factor that rises with distance, which is
impossible: material already deposited does not come back. This package gives
**0.0188** there, which is monotone and one decimal point from the printed 0.19.
The test leaves that cell out and pins the computed value instead.

That benchmark also fixes what a lid does to a source standing at it. Table 3.7
tabulates a 100 m release in category F, whose mixing depth is 100 m, and only
the capped form reproduces it — 0.356 and 0.089 at 50 and 100 km against a
published 0.33 and 0.083, where the uncapped form gives 0.599 and 0.311. What
happens to a plume *above* the lid is a separate matter, and one of convention;
see the mixing layer below.

Reproducing published benchmarks at all needed one thing the package lacked: they
state the effective release height and the transport wind speed as inputs, where
the package derives both. `PrescribedPlume` wraps a `Site` and takes them as
given — the two σ as well, for Turner's problems — bypassing plume rise, the wind
profile and the dispersion curves. A prescribed σ_z holds at one distance only,
so the depletion integral refuses it.

### The mixing layer

A plume does not disperse upwards forever. Turbulent mixing is capped by an
inversion, and once σ_z approaches that depth the plume is trapped between the
ground and the lid and reflects off both. The package had no lid, which is what
limited agreement with published depletion tables in the far field.

It has one now: `HPA-RPD-058` §3.2.2.1 Eq. (3.4), the image sum over virtual
sources at `2sA ± h_e`. The report truncates it at `|s| = 1` and changes to the
uniform profile `1/A` of Eq. (3.5) once σ_z reaches the depth, noting that the
series "can be summed to any prescribed accuracy". The package sums it: the
image series while `Σ_z < 0.7 A`, and above that the cosine series of the same
function, whose leading term *is* Eq. (3.5). There is no change of formula, the
activity between the ground and the lid integrates to one to within 4 × 10⁻¹⁶
at every σ_z, and the profile tends to `1/A` of its own accord: 1.4 % above it
at `Σ_z = A`, 5 × 10⁻⁹ at `2A`. The truncated sum stays within 0.06 % of the
converged one up to `Σ_z = A` and is 3 % low at `1.5 A`.

Depths are Table 3.5(a), which the report attributes to [Clarke1979](@citet)
and [Jones1980](@citet):

| A | B | C | D | E | F |
|---|---|---|---|---|---|
| 1300 m | 900 m | 850 m | 800 m | 400 m | 100 m |

Deeper in unstable air, shallower in stable, which is the physical ordering. The
depth is a site measurement where there is one, so `MixingLayer` takes any six
depths. In a configuration, `[mixing_layer] scheme` selects `tabulated` (the
default), `uniform` with a `uniform_depth` — HPA-RPD-058 §3.2.2.2.3 recommends
800 m "for all conditions" when the real one is unknown — `custom` with six
`depths`, or `unbounded`.

**What it changes**, as a ratio to the unbounded field at ground level, for the
reference case:

| | 10 km | 20 km | 50 km |
|---|---|---|---|
| A | 1.06 | **1.45** | **2.23** |
| B | 1.00 | 1.12 | 1.67 |
| C | 1.00 | 1.01 | 1.16 |
| D | 1.00 | 1.00 | 1.00 |
| E | 1.00 | 1.00 | 1.01 |
| F | **2.28** | **2.14** | **2.12** |

The unstable classes reach the lid first, despite having the deepest one,
because σ_z grows fastest there. The lid can only raise ground-level
concentration — it reflects back down material that would otherwise have kept
rising — so the correction is in the conservative direction.

**A plume above the lid.** Class F is a different case: its tabulated depth is
100 m while plume rise puts the effective release at 103.5 m. Published practice
treats a plume above the lid in two ways, and the package carries both, as
`above_lid`:

- `rise_inhibited`, the default. NRPB-R157 [Jones1983](@cite) §B2.3: "plume rise will be inhibited
  by a capping inversion to the mixing layer. If a plume rises into such an
  inversion the amount of material in the mixing layer, and hence ground-level
  concentration, will be reduced." Holding the whole plume at the lid,
  `H = min(H, A)`, is therefore the conservative reading. It is continuous in
  `H`, and it is the case Table 3.7 tabulates. It is also what doubles class F
  in the table above: a source at a reflecting lid is its own image.
- `full_penetration`. EPA ISC3, User's Guide vol. II [EPA1995](@cite) §1.1.6.1: "if the effective
  stack height, he, exceeds the mixing height, zi, the plume is assumed to fully
  penetrate the elevated inversion and the ground-level concentration is set
  equal to zero". ISC3 also takes stable air as unbounded, which is a depth of
  `inf` for classes E and F.

An earlier version of this package ignored the lid where `H > A`, on the
argument that such a plume is decoupled from the ground until the inversion
breaks — fumigation, which is not modelled here. That matched neither published
convention and made the field discontinuous: a plume 0.1 m below a 100 m lid
gave twice the ground-level concentration of one 0.1 m above it.

In the reference case the choice moves the undepleted long-term χ/Q at 20 km in
the most exposed sector from 4.2 × 10⁻⁹ s m⁻³ with no lid to 5.2 × 10⁻⁹ with
the rise inhibited, and to 3.4 × 10⁻⁹ with full penetration.

Note that the analytic invariants — the crosswind integral, the sector
average, the `Σ_z = H/√2` maximum — are properties of the *unbounded* Gaussian
and are asserted against a site with the lid switched off. They are statements
about the kernel's normalisation, not about the atmosphere.

### The depletion and resuspension parameters

These were the last constants in the package with no attribution. They have one
now.

**Resuspension.** `K = A exp(−λ₁t) + B exp(−λ₂t)` with
`A = 10⁻⁵ m⁻¹, B = 10⁻⁹ m⁻¹, λ₁ = 10⁻² d⁻¹, λ₂ = 2 × 10⁻⁵ d⁻¹` is **IAEA Safety
Series No. 57** [IAEA1982](@cite), §3.6, Eq. (3.14A) — all four constants, exactly, in the
same units, read in the primary document, which states them in that sentence.
It is reference [3] of NSR-23's own bibliography, which is how they got here.
Safety Series 57 is now superseded and says so on every page, but no successor
restates these constants. That report also brackets them: A over 10⁻⁶–10⁻⁴ and B over
10⁻¹⁰–10⁻⁸ m⁻¹, with the fast half-life "of the order of weeks" (this gives
69.3 d) and the slow one "in the range 50 to 100 years" (94.9 yr). Asserted.

A later model supersedes it — [MaxwellAnspaugh2011](@citet), also in
NUREG/CR-7270 [Bixler2022](@cite) — which keeps the amplitudes but puts the fast
decay constant at 0.07 d⁻¹ rather than 0.01. The package therefore runs high
with elapsed time: 1.8× at ten days, 6× at thirty. Recorded, not changed.

**Washout.** The table is **NSR-23 Table 7**, which the normative attributes to
CAN/CSA-N288.2-M91. It has **two species rows** and this package originally
carried only one; both are here now, selected by `WashoutSpecies`:

| at 1 mm/h | rain Λ_L | rain Λ_H | snow Λ_L | snow Λ_H |
|---|---|---|---|---|
| Tritium and iodine | 1 × 10⁻⁵ | 2 × 10⁻⁴ | 1 × 10⁻⁷ | 4 × 10⁻⁷ |
| All other nuclides | 2 × 10⁻⁵ | 3 × 10⁻⁴ | **5 × 10⁻⁴** | **2 × 10⁻²** |

The intensity dependence is the published one: fitting each column in log–log
gives an exponent of **0.753** for three of four, and `Λ ∝ J^0.75` is
[Slinn1977](@citet) via **NRPB-R322** [ADMLC2001](@cite) §3.1.1. Safety Series 57 §3.4.2 takes a
third position, `Λ = aI`, linear; NRPB-R157 §D3.4 brackets the exponent at 0.5
to 1.0, which contains both, so neither is asserted against the other.

The rain amplitudes bracket every published value found. At 1 mm/h the German
**AVV** Anhang 7 Tabelle 3 gives 7 × 10⁻⁵ s⁻¹ for aerosols and elemental iodine
and 3.5 × 10⁻⁵ for tritiated water; **Safety Series 57 Table II** gives
`a = 1.6 × 10⁻⁴` h(mm·s)⁻¹ for particulates and `1.1 × 10⁻⁴` for elemental
iodine, so 1.6 × 10⁻⁴ and 1.1 × 10⁻⁴ s⁻¹ at that rate. All four fall inside
**both** species rows.

**Snow is where it breaks down, and the normative disagrees with itself.** The
tritium row's snow values are the rain values divided by 100 and 500 — a
particle-scavenging suppression, which IAEA TECDOC-379 [IAEA1986](@cite) §3.5.4 supports for
*particles* (inorganic iodine in powder snow at 0.2 mm/h, 5 × 10⁻⁸ s⁻¹, against
1.7 × 10⁻⁵ for the same species in rain). But the other-nuclides row runs the
*other* way, three to four orders of magnitude **above** its own rain values.

Snow scavenging of tritiated water is isotopic exchange at the crystal surface,
not impaction, so the suppression is the wrong mechanism for it.
[Ogram1985](@citet) measured it, and the measurement lands almost exactly on
NSR-23's *other-nuclides* snow lower limit:

| mm/h | Ogram HTO | NSR-23 other-nuclides snow Λ_L | ratio |
|---|---|---|---|
| 0.5 | 2.88 × 10⁻⁴ | 3 × 10⁻⁴ | 0.96 |
| 1 | 4.20 × 10⁻⁴ | 5 × 10⁻⁴ | 0.84 |
| 3 | 7.78 × 10⁻⁴ | 8 × 10⁻⁴ | 0.97 |
| 5 | 1.04 × 10⁻³ | 1 × 10⁻³ | 1.04 |

Four rates, agreement within 16 %, from a measurement the normative does not
cite. That is the evidence that the tritium row's snow column is the anomaly.
`[model] washout = "hto"` selects Ogram's correlation; rain is unaffected.

### The tritium deposition velocities

The last parameter in the reference configuration without a citation, and it
turned out to have one. The notes to **NSR-23 Table 6** state the deposition
velocity of HTO as **0.4–0.8 × 10⁻² m/s** and that of HT as an order of magnitude
lower, **0.04–0.05 × 10⁻² m/s**, both attributed to experimental measurement —
[Murphy1993](@citet) — under the stated condition that the tropopause be taken at
12–15 km. Those are exactly the four numbers the package carries.

The lower bound is close to what others measure: 0.5 cm/s is the MACCS2 default
and 0.42 the Savannah River Site value, both inside the range. AECL's 0.392–0.444
straddles it — its lower end sits 2 % below the norm's 0.4 — so the normative's
lower bound is at the edge of the measurements rather than inside them, which is
asserted as such.

**That HTO deposits at all is a divergence, and a deliberate one.** IAEA SRS-19
§3.9, EUR 15760 [Simmonds1995](@cite) §3.2 and HPA-RPD-058 §3.2.2.3 each assign tritium a deposition
velocity of **zero**, handling it by specific activity instead. This package
follows the normative the work was done under and says so; setting both
velocities to zero follows the others.

### Where schemes disagree, the choice is in the configuration

Three plume-rise constants differ between published schemes, and the 2021 code
took a position on each without recording one. They are now a `RiseCoefficients`
value selected by `[model] plume_rise`, defaulting to Briggs:

| | Briggs, the default | 2021 code | Also published |
|---|---|---|---|
| Combined-law buoyancy denominator | **0.72** = 2β², β = 0.6 | 0.5 | — |
| Neutral momentum rise, `c w₀D/u` | **3** — [Briggs1969](@citet) Eq. 5.2, EPA ISC3 Eq. (1-16) | 1.5, the momentum term of [Holland1953](@citet), [Turner1970](@citet) Eq. (4.1) | — |
| Stable final rise, `c[F/(uS)]^(1/3)` | **2.6** — Briggs, in the Handbook on Atmospheric Diffusion [Hanna1982](@cite) | 2.6 | 2.4, NRC XOQDOQ |

The first is the one that matters. The combined momentum-and-buoyancy law must
reduce to the pure-buoyancy law when the momentum flux vanishes, and with 0.5 it
did not: it left `(6F x²/u³)^(1/3)` against the two-thirds law's
`(1.6³F x²/u³)^(1/3)`, an overshoot of **13.6 %**. The momentum half of the same
expression was already exactly Briggs — `β_j = 1/3 + u/w₀` — which placed the
discrepancy in the buoyancy half alone, and two independent routes give the
correction: `2β² = 0.72` and `3/1.6³ = 0.7324`.

With 0.72 the law reduces as it should, to within **0.6 %** — the residual being
that `(3/0.72)^(1/3) = 1.60915` where the literature rounds the two-thirds
coefficient to 1.6, the same rounding that makes the neutral final rise 21.425
against a published 21.4.

`plume_rise = "nsr23"` restores the old constants exactly, and the
fidelity tests use it, so the 2021 results remain reproducible.

**The reference case does not move**, and that is worth saying rather than
leaving as a surprise: the CANDU stack is buoyancy-dominated — 57.7 m of buoyant
rise against 10.6 m of momentum rise in class D, an imbalance far outside the
tolerance for the combined law — so it takes the pure-buoyancy branch, which
none of these constants touch. The correction bites on momentum-dominated and
balanced releases, which is where the 13.6 % was.

Selecting `xoqdoq` does move it: the reference atmosphere is stably stratified,
so the stable branch is active, and the effective release height at 2 km in
class D falls from 108.0 m to 103.6 m.

Resuspension is selectable the same way: `[model] resuspension` takes
`iaea_ss57` (the default, Safety Series 57) or `maxwell_anspaugh` (2011, also
NUREG/CR-7270).

**The building-wake coefficient now defaults to 1.0**, which is IAEA SRS-19
Eq. (6), `Σ_z = (σ_z² + A_B/π)^(1/2)`, and the German AVV Eqs. (4.31)/(4.32),
`√(σ² + I_G²/π)`. NSR-23 used 1.5 and no source outside it was found;
`NSR23_WAKE_COEFFICIENT` restores it. RG 1.111 Eq. (9) is a third convention
again, applying 0.5 to the building *height* rather than its area.

**Snow washout is now selectable by species.** The normative's snow columns are
the rain columns divided by exactly 100 and 500 — a particle-scavenging
suppression, and the right physics for aerosols and reactive gases: IAEA
TECDOC-379 §3.5.4 gives 5 × 10⁻⁸ s⁻¹ for inorganic iodine in powder snow at
0.2 mm/h against 1.7 × 10⁻⁵ for the same species in rain, a factor of some 340
the same way. But snow scavenging of **tritiated water is isotopic exchange at
the crystal surface, not impaction**, and Ogram measures it about **a thousand times above** the normative column, not below it.
`[model] washout = "hto"` selects Ogram's correlation for snow; rain is
identical either way. The reference case of this package is HTO, so the choice
matters for it.

The third row is worth a note. This code writes the neutral final rise as
`1.6F^(1/3)(3.5x_f)^(2/3)/u`, which does not look like the published
`21.4F^(3/4)/u`. Substituting `x_f = 14F^(5/8)` collapses it: `1.6·49^(2/3) =
21.425`, and `1.6·119^(2/3) = 38.71` for the other branch. They are one
expression, and the literature rounds the constant. The tests assert both the
algebra and the numbers.

σ_z is not Briggs — it is Hosker, as above, carrying an explicit roughness
correction that Briggs does not have. So the package mixes **Briggs σ_y with
Hosker σ_z**, a pairing no single publication uses, and the two are therefore
asserted against different sources. Against Briggs, σ_z is only bracketed: the
ratio runs from **0.41** (class B at 10 km) to **1.49** (class E at 10 km) over
0.1–10 km, a factor of 2.4 taken symmetrically, and it is systematically ordered
— less vertical spread than Briggs in unstable air, more in stable, closest in
neutral. That is the spread that separates published σ schemes from one another,
and `figures/validation/dispersion_parameters.png` is the picture of it.

The last check needs a word. `Σ_z = H/√2` is derived holding `H` fixed and
`Σ_y ∝ Σ_z`; asserted under those assumptions it is exact. In the real field the
maximum lands 5 % away from it, and the cause is not error but `Σ_y/Σ_z`
drifting from 1.94 to 2.43 across the peak region while plume rise has already
saturated. All six Pasquill classes and three sector counts are covered.
