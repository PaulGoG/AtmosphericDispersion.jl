#=
Tabulated coefficients of the dispersion parameterisation.

These are fixed tables of the governing normative — CNCAN NSR-23, the Romanian
nuclear regulator's norm the 2021 thesis was written against — reproduced as
Tables 1-4 and 7 of the thesis. They are compile-time constants rather than data files: they
are small, they never vary between runs, and the values are read inside the
innermost loops of every field evaluation, where a keyed lookup into a data
frame costs more than the dispersion calculation itself.

The vertical dispersion tables are Hosker's, and the normative reproduces them:
Hosker, R.P. Jr., "Estimates of dry deposition and plume depletion over forests
and grassland", IAEA-SM-181/19, pp. 291-309 in Physical Behaviour of Radioactive
Contaminants in the Atmosphere, IAEA, Vienna (1974) — an analytic fit to F.B.
Smith (1972) and Briggs (1973). They are printed in Smith and Simmonds (eds.),
HPA-RPD-058, Health Protection Agency (2009), Table 3.3, and in Clarke, R.H.,
NRPB-R91 (1979), Table 3, and both printings are asserted in the test suite.

Indexing is positional, through `classindex` and `Int(::enum)`, so every lookup
is a tuple index with a concrete return type.
=#

"""
    VerticalShapeCoefficients

Coefficients of the shape function of the vertical dispersion parameter,

    g(x) = a₁ x^b₁ / (1 + a₂ x^b₂)

with `x` in metres. Keyed by Pasquill class; see [`vertical_shape_coefficients`](@ref).
"""
struct VerticalShapeCoefficients
    a₁::Float64
    b₁::Float64
    a₂::Float64
    b₂::Float64
end

const _VERTICAL_SHAPE = (
    VerticalShapeCoefficients(0.112, 1.06, 5.38e-4, 0.815),    # A
    VerticalShapeCoefficients(0.130, 0.95, 6.52e-4, 0.750),    # B
    VerticalShapeCoefficients(0.112, 0.92, 9.05e-4, 0.718),    # C
    VerticalShapeCoefficients(0.098, 0.889, 1.35e-3, 0.688),   # D
    VerticalShapeCoefficients(0.0609, 0.895, 1.96e-3, 0.684),  # E
    VerticalShapeCoefficients(0.0638, 0.783, 1.36e-3, 0.672),  # F
)

"""
    vertical_shape_coefficients(class)

[`VerticalShapeCoefficients`](@ref) of the given Pasquill `class`.
"""
vertical_shape_coefficients(class::PasquillClass) = _VERTICAL_SHAPE[classindex(class)]

"""
    RoughnessCoefficients

Roughness length `z₀` in metres and the coefficients of the roughness
correction factor applied to the vertical dispersion parameter. See
[`roughness_coefficients`](@ref) and [`roughness_correction`](@ref).
"""
struct RoughnessCoefficients
    z₀::Float64
    c₁::Float64
    d₁::Float64
    c₂::Float64
    d₂::Float64
end

# HPA-RPD-058 Table 3.3, second panel; identical in NRPB-R91 Table 3.
#
# NSR-23 Table 2 prints 1.58 and 2.08 for the first two, which is what the 2021
# code faithfully carried. The error is the normative's, not the thesis's: it
# attributes its Table 2 to CAN/CSA-N288.2-M91 and UNSCEAR 2000, and somewhere
# in that chain Hosker's 1.56 and 2.02 became 1.58 and 2.08. Following the
# published values makes σ_z smaller by 1.9 % over grassland and water and 3.6 %
# over arable land.
const _ROUGHNESS = (
    RoughnessCoefficients(0.01, 1.56, 0.048, 6.25e-4, 0.45),    # grassland and water
    RoughnessCoefficients(0.04, 2.02, 0.0269, 7.76e-4, 0.37),   # arable
    RoughnessCoefficients(0.10, 2.72, 0.0, 0.0, 0.0),           # pasture
    RoughnessCoefficients(0.40, 5.16, -0.098, 18.6, -0.225),    # rural
    RoughnessCoefficients(1.00, 7.37, -0.0957, 4.29e3, -0.6),   # forest and urban
    RoughnessCoefficients(4.00, 11.7, -0.128, 4.59e4, -0.78),   # metropolis
)

"""
    roughness_coefficients(roughness)

[`RoughnessCoefficients`](@ref) of the given [`RoughnessClass`](@ref).
"""
roughness_coefficients(roughness::RoughnessClass) = _ROUGHNESS[Int(roughness)]

"""
    roughness_length(roughness)

Roughness length `z₀` of the given class, in metres.
"""
roughness_length(roughness::RoughnessClass) = roughness_coefficients(roughness).z₀

# Table 3: coefficient of the lateral dispersion parameter, by Pasquill class.
const _LATERAL_COEFFICIENT = (0.22, 0.16, 0.11, 0.08, 0.06, 0.04)

"""
    lateral_coefficient(class)

Coefficient `c₃` of the lateral dispersion parameter for the given Pasquill
class; see [`lateral_dispersion`](@ref).
"""
lateral_coefficient(class::PasquillClass) = _LATERAL_COEFFICIENT[classindex(class)]

# Table 4: exponent of the power-law wind profile, indexed [surface, class].
const _PROFILE_EXPONENT = (
    (0.03, 0.05, 0.06, 0.08, 0.10, 0.12),   # water
    (0.10, 0.15, 0.20, 0.25, 0.35, 0.40),   # agricultural
    (0.16, 0.24, 0.32, 0.40, 0.56, 0.64),   # forest and urban
)

"""
    profile_exponent(surface, class)

Exponent `m` of the power-law wind profile `u(z) = u₁₀ (z/10)^m` for the given
[`WindProfileSurface`](@ref) and Pasquill class.

The exponent grows with both surface roughness and atmospheric stability: shear
is strongest over rough ground under a stable stratification.
"""
profile_exponent(surface::WindProfileSurface, class::PasquillClass) =
    _PROFILE_EXPONENT[Int(surface)][classindex(class)]

"""
    WashoutCoefficients

Lower and upper washout coefficients `Λ` in s⁻¹ for a given precipitation type
and intensity. The pair brackets the reported range; `low` is used where an
underestimate of removal is conservative for the airborne concentration and
`high` where an overestimate of deposition is.
"""
struct WashoutCoefficients
    low::Float64
    high::Float64
end

"""
    PRECIPITATION_RATES

The precipitation intensities in mm/h at which the washout coefficients are
tabulated. [`washout_coefficients`](@ref) is defined at these values and
nowhere else.
"""
const PRECIPITATION_RATES = (0.5, 1.0, 3.0, 5.0)

"""
    WashoutSpecies

Which row of the normative's washout table applies. NSR-23 Table 7, which it
attributes to CAN/CSA-N288.2-M91, tabulates two:

  - `WASHOUT_TRITIUM_IODINE` — tritium and iodine, the default, and the case
    this package's reference configuration is.
  - `WASHOUT_OTHER_NUCLIDES` — every other radionuclide. Its rain values are
    roughly twice the tritium row's; its **snow** values are three to four
    orders of magnitude larger, not smaller.

That the two snow rows differ by so much in opposite directions is the reason
the snow scavenging of tritiated water is worth stating explicitly rather than
inheriting. See [`WashoutModel`](@ref).
"""
@enum WashoutSpecies begin
    WASHOUT_TRITIUM_IODINE = 1
    WASHOUT_OTHER_NUCLIDES = 2
end

# NSR-23 Table 7, both rows, at PRECIPITATION_RATES. Indexed
# [species][precipitation][rate].
const _WASHOUT = (
    (   # tritium and iodine
        (   # rain
            WashoutCoefficients(5.0e-6, 1.0e-4),
            WashoutCoefficients(1.0e-5, 2.0e-4),
            WashoutCoefficients(2.0e-5, 4.0e-4),
            WashoutCoefficients(3.0e-5, 6.0e-4),
        ),
        (   # snow
            WashoutCoefficients(5.0e-8, 2.0e-7),
            WashoutCoefficients(1.0e-7, 4.0e-7),
            WashoutCoefficients(2.0e-7, 8.0e-7),
            WashoutCoefficients(3.0e-7, 1.0e-6),
        ),
    ),
    (   # all other radionuclides
        (   # rain
            WashoutCoefficients(1.0e-5, 2.0e-4),
            WashoutCoefficients(2.0e-5, 3.0e-4),
            WashoutCoefficients(3.0e-5, 7.0e-4),
            WashoutCoefficients(5.0e-5, 1.0e-3),
        ),
        (   # snow
            WashoutCoefficients(3.0e-4, 1.0e-2),
            WashoutCoefficients(5.0e-4, 2.0e-2),
            WashoutCoefficients(8.0e-4, 4.0e-2),
            WashoutCoefficients(1.0e-3, 5.0e-2),
        ),
    ),
)

"""
    WashoutModel

Which scavenging scheme to use.

  - `WASHOUT_NORMATIVE` — the tabulated scheme of the normative, and the
    default. Its snow columns are the rain columns divided by exactly 100 and
    500, a suppression appropriate to **particles and reactive gases**: IAEA
    TECDOC-379 §3.5.4 gives 5 × 10⁻⁸ s⁻¹ for inorganic iodine in powder snow at
    0.2 mm/h against 1.7 × 10⁻⁵ for the same species in rain, a factor of some
    340 in the same direction.
  - `WASHOUT_HTO` — the same rain columns, with snow from Ogram, *Precipitation
    Scavenging of Tritiated Water Vapour (HTO)*, Ontario Hydro Research Division
    85-233-K (1985), §6.0 Eq. (38):

        Λ_s = 1.2×10⁻⁴ R^0.33 + 3.0×10⁻⁴ R^0.64   s⁻¹,  R in mm/h

    Snow scavenging of tritiated water is **isotopic exchange at the crystal
    surface**, not impaction, so the particle suppression is the wrong physics
    for it. Ogram measures scavenging about three orders of magnitude *above* the
    normative snow column, and calls the correlation an approximate upper limit.

The reference case of this package is HTO, so the two differ for it by a factor
of roughly a thousand in snow. They are identical in rain.
"""
@enum WashoutModel begin
    WASHOUT_NORMATIVE = 1
    WASHOUT_HTO = 2
end

"""
    ogram_snow_washout(rate)

Washout coefficient in s⁻¹ for tritiated water in snow at `rate` mm/h of water
equivalent, from Ogram (1985) Eq. (38).

The report's own Table V reproduces this at 0.5, 1 and 2 mm/h but prints
2.0 × 10⁻⁴ at 0.1 mm/h where the equation gives 1.25 × 10⁻⁴. That inconsistency
is in the original; the equation is used here.
"""
function ogram_snow_washout(rate::Real)
    rate ≥ 0 || throw(DomainError(rate, "precipitation rate cannot be negative"))
    return 1.2e-4 * rate^0.33 + 3.0e-4 * rate^0.64
end

"""
    washout_coefficients(precipitation, rate, model = WASHOUT_NORMATIVE)

[`WashoutCoefficients`](@ref) for the given precipitation type and intensity in
mm/h, under the chosen [`WashoutModel`](@ref).

The table is defined only at the intensities in [`PRECIPITATION_RATES`](@ref);
any other value throws. The tabulated points follow a power law in intensity
closely enough that interpolating between them would be defensible, but that is
a modelling decision rather than a lookup, and it is not made here.
"""
function washout_coefficients(
    precipitation::PrecipitationType,
    rate::Real,
    model::WashoutModel = WASHOUT_NORMATIVE,
    species::WashoutSpecies = WASHOUT_TRITIUM_IODINE,
)
    i = findfirst(==(float(rate)), PRECIPITATION_RATES)
    isnothing(i) && throw(
        ArgumentError(
            "the washout table is defined at intensities $(PRECIPITATION_RATES) mm/h, got $rate",
        ),
    )
    tabulated = _WASHOUT[Int(species)][Int(precipitation)][i]
    if model == WASHOUT_HTO && precipitation == PRECIPITATION_SNOW
        # Ogram gives one correlation, described as an upper limit, so it is
        # taken as the high bound and the tabulated value as the low one.
        return WashoutCoefficients(tabulated.low, ogram_snow_washout(rate))
    end
    return tabulated
end
