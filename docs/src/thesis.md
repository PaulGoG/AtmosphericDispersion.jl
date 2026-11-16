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
| Plume-rise constants | none: the stack is buoyancy-dominated, 55.7 m of buoyant rise against 10.6 m of momentum, so the corrected constants are never reached |
| Building-wake coefficient | none: no buildings in the reference configuration |
| Roughness coefficients | none: the corrected rows are z₀ = 0.01 and 0.04 m, and the case is pasture at 0.1 m |
| Mixing layer, which the 2021 code did not have | none within 2 km; in the most exposed sector +2 % at 5 km, +11 % at 10 km and +22 % at 20 km, nearly all of it class F, whose plume is held at its 100 m lid |
| Stable final rise on the mean wind over the rise | effective heights lower by 1 to 3 m; χ/Q in the most exposed sector higher by 5 % at 1 km, 4 % at 2 km and 0.3 % at 20 km |

The published dose maps are rotated by a half turn, and the most-exposed sector
carries 1.51 times what was reported, more beyond a few kilometres, where the
lid adds up to 22 %. The remaining corrections affect only other configurations
(a momentum-dominated stack, a site with buildings, grassland or arable
roughness).

Both the 2021 code and this package apply the stable limit of the plume rise in
every stability class of this case; the [known limitations](limitations.md)
page gives the size of that effect, which is larger than any correction above
except the convention.
