# The 2021 thesis code

The [`original`](https://github.com/PaulGoG/AtmosphericDispersion.jl/tree/original) branch holds the code as submitted
for the BSc thesis at the Faculty of Physics, University of Bucharest, in June
2021, its files under English names since 2026, together with the thesis itself
(`BSc_thesis_2021.pdf`, in Romanian). `main` shares no history with it.

## How far the 2021 results move

Only some of the corrections recorded in the [validation](validation.md) page
and in the changelog affect the CANDU case the thesis assessed. Running that
configuration both ways:

| Correction | Effect on the reference case |
|---|---|
| Wind-direction convention | Most-exposed sector S instead of N (a half turn) and χ/Q there higher by a factor of 1.51 at every distance |
| Stable limits of the plume rise in the stable classes only, where the 2021 code applied them in every class | effective heights of classes A to D 130 to 153 m instead of 106 to 112 m; χ/Q in the most exposed sector lower by a factor of 0.48 at 1 km, 0.64 at 2 km, 0.89 at 5 km and 0.99 at 20 km |
| Plume-rise constants | none: the stack is buoyancy-dominated, 80.2 m of buoyant rise against 11.7 m of momentum in class D, so the corrected constants are never reached |
| Building-wake coefficient and rules | none: no buildings in the reference configuration |
| Roughness coefficients | none: the corrected rows are z₀ = 0.01 and 0.04 m, and the case is pasture at 0.1 m |
| Mixing layer, which the 2021 code did not have | none within 2 km; in the most exposed sector +2 % at 5 km, +11 % at 10 km and +23 % at 20 km, nearly all of it class F, whose plume is held at its 100 m lid |
| Stable final rise on the mean wind over the rise | effective heights of classes E and F lower by 2.6 and 2.7 m; χ/Q in the most exposed sector higher by 1 % at 2 and 5 km and 0.2 % at 20 km |
| Washout and wet deposition | none: the reference configuration has no precipitation |

The published dose maps are rotated by a half turn, which moves their maximum
without changing its value. Peak against peak, the corrected maximum is about
half the published one at 1 km, two thirds at 2 km, 0.92 of it at 5 km and
1.2 times it at 20 km. Near the source the 2021 code held the plume of the
unstable and neutral classes down with the stable limit of the rise; further
out the lid adds up to 23 %. The remaining corrections affect only other
configurations (a momentum-dominated stack, a site with buildings, grassland or
arable roughness, precipitation).
