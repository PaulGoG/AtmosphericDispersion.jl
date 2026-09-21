# The 2021 thesis code

The [`original`](https://github.com/PaulGoG/AtmosphericDispersion.jl/tree/original) branch holds the code exactly as submitted
for the BSc thesis at the Faculty of Physics, University of Bucharest, in June
2021, together with the thesis itself (`BSc_thesis_2021.pdf`, in Romanian).
`main` shares no history with it.

## How far the 2021 results move

Every correction above is real, but only some of them bite on the CANDU case the
thesis assessed. Running that configuration both ways:

| Correction | Effect on the reference case |
|---|---|
| **Wind-direction convention** | **Most-exposed sector S instead of N — a half turn — and χ/Q there higher by a factor of 1.51 at every distance** |
| Plume-rise constants | none: the stack is buoyancy-dominated, 57.7 m of buoyant rise against 10.6 m of momentum, so the corrected constants are never reached |
| Building-wake coefficient | none: no buildings in the reference configuration |
| Roughness coefficients | none: the corrected rows are z₀ = 0.01 and 0.04 m, and the case is pasture at 0.1 m |
| Mixing layer | none within 30 km in class D; σ_z never approaches 800 m there |

So the answer is a single number: **the published dose maps are rotated by a half
turn, and the most-exposed sector carries 1.51 times what was reported.**
Everything else corrected here would change a different configuration — a
momentum-dominated stack, a site with buildings, grassland or arable roughness,
or the far field in unstable air — and leaves this one alone.

That is worth saying plainly because the reverse would have been easy to assume.
A dozen corrections do not compound into a dozen shifts; most of them are
inactive in any one configuration, and the one that dominated was the one that
needed no numerics to see.
