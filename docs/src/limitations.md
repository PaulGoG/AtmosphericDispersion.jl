# Known limitations

An independent review of the physics against the literature, in November 2026,
found the points below. They are not corrected in the 1.0 series, because each
correction changes results, and several of the rules come from the normative
the 2021 work followed. The numbers refer to the reference configuration.

## One stability parameter for all classes

`Site` forms the stability parameter `S` from the single `lapse_rate` of its
`Atmosphere` and applies it in all six Pasquill classes. With a positive `S`
the stable limits of the plume rise therefore act in the unstable and neutral
classes as well. EPA ISC3 [EPA1995](@cite) §1.1.4.6 and XOQDOQ
[Sagendorf1982](@cite) compute a stability parameter for the stable classes
only and use the unstable–neutral rise in classes A to D.

The reference lapse rate, 0.02 K/m, is a strong inversion, and the stable limit
is the smallest in every class:

| Class | Stable limit [m] | Neutral limit [m] | Effective height [m] | With the neutral limit [m] |
|---|---|---|---|---|
| A | 61.6 | 102.2 | 111.9 | 152.5 |
| B | 59.5 | 94.3 | 109.8 | 144.6 |
| C | 57.6 | 87.0 | 107.9 | 137.3 |
| D | 55.7 | 80.2 | 106.0 | 130.5 |
| E | 52.1 | 68.2 | 102.4 | |
| F | 50.5 | 63.0 | 100.8 | |

The lower effective heights of classes A to D raise the ground-level
concentration. Against a calculation that applies the stable limits in classes
E and F only, the long-term χ/Q of the most exposed sector is higher by a
factor of 2.1 at 1 km, 1.6 at 2 km and 1.13 at 5 km, and by 1.4 % at 20 km. The
reference results are therefore conservative in the near field. A run that
needs the unstable–neutral rise in every class can set `lapse_rate` to the dry
adiabatic value, −g/c_p, which makes `S` zero; no single value is consistent
with all six classes.

## Washout

The wet depletion factor is `exp(−Λ t_w)` at every distance, with `t_w` the
duration of the precipitation event. Material that reaches the distance `x` has
been in the rain for at most the travel time `x/u`, which is the time
HPA-RPD-058 [SmithSimmonds2009](@cite) Eq. (3.13) uses. For an hour of rain at
1 mm/h the factor is 0.965 everywhere, where the travel time gives at most
0.998 at 1 km and 0.967 at 20 km, so the air concentration near the source is
underestimated. The same factor enters the long-term average as if the event
lasted the whole release.

Wet deposition is formed from the whole released activity whatever the duration
of the event. It does not vanish for an event of zero duration, and its
integral over the ground is not bounded by the release: within 20 km it is 0.64
of the released activity at 1 mm/h and 1.8 times it at 5 mm/h.

## Bracketing values and the activity budget

Depletion takes the low deposition velocity and the low washout coefficient,
and deposition takes the high ones. Each choice is conservative for its own
pathway, but together they do not conserve activity: the ground receives more
than the plume loses, by the ratio of the two values. For the dry deposition of
HTO that ratio is 2.

## Dry deposition of tritiated water

A deposition velocity above zero depletes the plume and lowers the air
concentration downwind. Deposited HTO is re-emitted, and
[Lee2012](@citet) conclude that an effective deposition velocity of zero is
the conservative choice for safety analyses, after the 0.5 cm/s of the MACCS2
guidance had been questioned as not bounding. IAEA SRS-19, EUR 15760 and
HPA-RPD-058 also take zero. The package follows the normative; setting both
velocities to zero gives the conservative air concentration.

## Buildings

The rules below are those of the normative as the 2021 code implemented them.

- The equivalent building is a mean over the buildings near the stack, weighted
  by inverse distance. A small structure close to the stack therefore lowers
  the equivalent height of a large one further away: a 3 m shed at 5 m reduces
  a 60 m building at 50 m to an equivalent 8.2 m, and the release is no longer
  entrained. IAEA SRS-19 §3.3 takes the building that most influences the flow.
- A release entrained into the cavity leaves from ground level but still
  receives the full plume rise, computed with the wind at 10 m. For a 60 m
  building 25 m from the reference stack the effective height in class D is
  59.3 m. SRS-19 §3.5 takes an effective height of zero in the wake.
- The downwash and entrainment tests use the wind at stack height averaged over
  the six classes, so every class has the same release height, and the
  entrainment test changes the release height discontinuously at its 5 m/s
  threshold.

## Plume rise

- Where the two final rises agree within 10 % the combined law applies, capped
  at their sum; otherwise the larger applies alone. The final rise is therefore
  discontinuous at the edge of that band, by nearly a factor of two.
- The transitional momentum law is written for an efflux of ambient density,
  while the combined law takes the density ratio through the momentum flux; for
  the reference density ratio of 0.5 the two differ by a factor of 1.26. The
  reference case is dominated by buoyancy and is not affected.
- `exit_temperature` is not used: the buoyancy flux is formed from the exit
  density alone, and the two are not checked against each other.

## Eimutis–Konicek at large spread

XOQDOQ and PAVAN limit both dispersion parameters of this scheme to 1000 m
[Sagendorf1982, Bander1982](@cite). The package evaluates the fit without the
limit. Its σ_z exceeds 1000 m beyond 1.5 km in class A and beyond 7.6 km in
class B, and its σ_y beyond 6.4 km in class A to 38 km in class F, all inside
the stated range of validity. Under the default mixing layer the plume already
fills the layer at those distances, and applying the limit would raise a
ground-level value of the sector-averaged forms by 10 % in class A and by less
than 1 % in the others. The short-release form divides by σ_y and changes by
the ratio of σ_y to 1000 m.

## Sector-averaged form off the axis

`dilution_extended` evaluates the plume at the distance along the axis of the
sector, and `dilution_long_term` at the radial distance. At equal radius, 10°
off the axis, the two differ by 0.2 % at 2 km and 2 % at 10 km.

## Low wind speed

The neutral limits of the plume rise vary as `1/u` and have no floor, so the
rise is unbounded as the wind speed tends to zero (4.8 km at 0.1 m/s for the
reference stack). EPA ISC3 does not let the wind at stack height fall below
1 m/s. The configuration loader requires a positive reference wind speed and
applies no further bound.
