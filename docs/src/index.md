```@meta
CurrentModule = AtmosphericDispersion
```

# AtmosphericDispersion

Gaussian-plume dispersion of stack releases: dilution factor, time-integrated
air concentration, dry and wet ground deposition, and resuspension, over
Pasquill–Gifford stability classes and a wind rose.

The solver core is nuclide-agnostic. A species enters only through its decay
constant and its deposition velocity, so tritium from a heavy-water reactor —
the case the original work addressed — is a parameterisation rather than an
assumption.

## A worked run

```julia
using AtmosphericDispersion

source = StackSource(;
    height = 50.3, diameter = 2.33, exit_velocity = 10.0,
    exit_density = 0.6, exit_temperature = 324.0,
)
atmosphere = Atmosphere(;
    reference_speed = 4.0, temperature = 287.0, density = 1.2,
    lapse_rate = 0.02, surface = SURFACE_AGRICULTURAL,
    roughness = ROUGHNESS_PASTURE,
)
site = Site(; source, atmosphere)

# A wind rose must declare which way its directions are read.
grid = SectorGrid(16)
rose = WindRose(grid, fill(1/16, 16), BlowingFrom(); stability = fill(1/6, 6))

# Long-term dilution 5 km due south of the stack, with tritium depletion.
dilution_long_term(0.0, -5000.0, site, rose; nuclide = TRITIATED_WATER)
```

Or from a configuration file:

```julia
config = load_configuration("config/reference.toml")
dilution_long_term(0.0, -config.extent, config.site, config.rose;
                   nuclide = config.nuclide)
```

## Conventions

Receptor positions are given in the geographic frame as displacements from the
stack, east first then north, in metres. Bearings are in radians clockwise from
north. A wind direction is the direction the wind blows **from**, as in
meteorological practice and in ADMS.

Sectors are indexed clockwise from north with sector 1 **centred** on north, so
that for sixteen sectors the index runs N, NNE, NE and so on. Centring rather
than starting at north puts the cardinal directions at sector centres, where
binning cannot be decided by floating-point round-off.

The long-term dilution formula weights by the frequency of wind blowing
**towards** the receptor's sector, which is the opposite reading. [`WindRose`](@ref)
therefore takes the convention as a required argument and stores
blowing-towards internally, so the formulae cannot be fed the wrong sense.
Getting this wrong rotates the whole field by half a turn while leaving every
magnitude and sum rule intact.

## The module

```@docs
AtmosphericDispersion
```

## Pages

```@contents
Pages = ["api/geometry.md", "api/parameters.md", "api/fields.md", "api/configuration.md"]
Depth = 2
```
