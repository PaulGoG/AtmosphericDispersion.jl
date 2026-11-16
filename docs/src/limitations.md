# Known limitations

An independent review of the physics against the literature, in November 2026,
and the published text of CNCAN NSR-23 [CNCAN2004](@cite) led to the
corrections of 1.1.0: the stable limits of the plume rise in the stable classes
only, washout for the travel time, a bounded wet deposition, the governing
building and the limit on wake broadening, the 1000 m limit of the
Eimutis–Konicek scheme, the cap of the combined rise and the minimum wind
speed. The [changelog](https://github.com/PaulGoG/AtmosphericDispersion.jl/blob/main/CHANGELOG.md)
lists them with their effect on the reference case. The points below remain.

## Bracketing values and the activity budget

Depletion of the airborne plume takes the low deposition velocity and the
low washout coefficient, and deposition the high ones, as NSR-23 Tables 6 and
7 direct. Each choice is conservative for its own pathway. Wet deposition
depletes its own column with the high coefficient, so the activity it deposits
cannot exceed the activity released in the precipitation. Dry deposition is
the high velocity times an air concentration depleted with the low one, so the
ground receives up to the ratio of the two velocities times what the plume
loses: 2 for HTO.

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
