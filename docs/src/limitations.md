# Known limitations

An independent review of the physics against the literature, in November 2026,
and the published text of CNCAN NSR-23 [CNCAN2004](@cite) led to the
corrections of 1.1.0: the stable limits of the plume rise in the stable classes
only, washout for the travel time, a wet deposition on the share of the release
in the precipitation, the governing building and the limit on the wake of a
ground-level release, the 1000 m limit of the Eimutis–Konicek scheme, the cap
of the combined rise and the minimum wind speed. The [changelog](https://github.com/PaulGoG/AtmosphericDispersion.jl/blob/main/CHANGELOG.md)
lists them with their effect on the reference case. The points below remain.

## Bracketing values and the activity budget

Depletion of the airborne plume takes the low deposition velocity and the
low washout coefficient, and deposition the high ones, as NSR-23 Tables 6 and
7 direct. Each choice is conservative for its own pathway. Dry deposition is
the high velocity times an air concentration depleted with the low velocity
and the low washout coefficient, so the ground receives up to the ratio of the
two velocities times what dry deposition takes from the plume, 2 for HTO. Wet
deposition is the largest deposit the precipitation can lay at each distance,
undepleted on the way, so its integral over the ground, `f A Λ_H X/u` out to
`X`, exceeds the activity released into the precipitation beyond
`X = u/Λ_H`, 10 km for HTO in rain at 5 mm/h and 6 m/s. The budgets are
separate: material that wet deposition removes is still in the air that feeds
dry deposition, so with precipitation the two deposits together can exceed
what the plume loses.

## Washout of a short precipitation

Washout of the airborne plume takes the material released into the
precipitation to stay in it for its whole travel time, the limit of a
precipitation much longer than that. For a shorter precipitation fixed over
the site the result depends on when it falls, which a configuration does not
state: for `Λ = 10⁻³ s⁻¹`, 900 s of precipitation in a 3600 s release and
20 km at 6 m/s, the surviving fraction lies between 0.47 and 0.92 by timing,
against 0.76 from the form. Where `Λ x/u` is large the precipitation must
outlast the travel time by far, `t_w ≫ (exp(Λ x/u) − 1)/Λ`.
`load_configuration` warns when a precipitation beginning with the release,
which leaves the most airborne, would leave more than 10 % more than the form
at the edge of the receptor grid. Setting no washout bounds the air
concentration for any timing. Wet deposition is not affected: it is the bound
for any timing.

## Dry deposition of tritiated water

A deposition velocity above zero depletes the plume and lowers the air
concentration downwind. Deposited HTO is re-emitted, and [Lee2012](@citet)
conclude that an effective deposition velocity of zero is the conservative
choice for safety analyses, after the 0.5 cm/s of the MACCS2 guidance had been
questioned as not bounding. IAEA SRS-19, EUR 15760 and HPA-RPD-058 also take
zero. The package follows NSR-23, whose notes to Table 6 give 0.4–0.8 cm/s;
setting both velocities to zero gives the conservative air concentration.

## Buildings

These are the rules of NSR-23 §6.2, which the package follows.

- A release below the top of the governing building is entrained into the
  cavity and released at ground level, and the plume rise is then added, as
  NSR-23 prescribes. For a 60 m building 25 m from the reference stack the
  effective height in class D is 120 m. The buoyancy flux is not lost when the
  plume is mixed into the cavity, but the rise of a plume that leaves a cavity
  is smaller than that of one leaving a stack; IAEA SRS-19 takes the release
  height as zero throughout the wake zone, which it calls pessimistic.
- Between one and two and a half building heights, the entrainment test
  depends on whether the mean wind at the stack, averaged over the six
  classes, exceeds 5 m/s, so the release height changes discontinuously at
  that wind speed.
- The `√3` limit on the wake broadening applies to a plume at effective
  height zero, a release at ground level with no plume rise, the case for
  which NRC XOQDOQ and NSR-23 §5 state it. A plume with any rise, including one
  entrained into the cavity and released at ground level, is broadened without
  limit, as NSR-23 Eq. (9) prescribes, since for it the limit would lower the
  ground-level concentration. The limit binds while `σ < √(C A/2π)`: for the
  60 m building of 3600 m² beside the reference stack, out to 0.54 km in class D
  and 2.6 km in class F. There a stack with the smallest efflux, whose plume
  rises, gives a lower concentration than one with none, at 200 m by a factor
  of 2.0 in class D and 5.1 in class F.
- The vicinity of a building is judged from the point given for it, which
  should be the point of the building nearest the stack; the package carries
  no footprint.

## Plume rise

- Where the two final rises agree within 10 % the combined law applies,
  otherwise the law of the larger mechanism alone. The far-field rise is
  continuous at the edge of that band; the near-field rise is not.
- The transitional momentum law is NSR-23 Eq. (21), written for an efflux of
  ambient density, while the combined law takes the density ratio through the
  momentum flux; for the reference density ratio of 0.5 the two differ by a
  factor of 1.26. The reference case is dominated by buoyancy.
- `exit_temperature` is not used: the buoyancy flux is formed from the two
  densities, as in NSR-23 Eq. (17.b), and the two are not checked against each
  other. The reference pair, 0.6 kg/m³ at 324 K, is lighter than steam at
  that temperature.
- The stable classes E and F share the one temperature gradient of the
  `Atmosphere`. EPA ISC3 takes 0.020 K/m of potential temperature for E and
  0.035 K/m for F where no measurement is available; the reference gradient,
  0.02 K/m of temperature, is 0.030 K/m of potential temperature.

## Meander

The broadening of `σ_y` for a release longer than ten minutes, `(t_R/600)^0.2`,
is defined in NSR-23 Eq. (8) between ten minutes and one hour, for which the
normative recommends the sector-averaged form beyond. The package applies it at
every duration. The Handbook on Atmospheric Diffusion gives an exponent of 0.25
to 0.3 between 1 h and 100 h.
