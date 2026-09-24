#=
Tabulated coefficients of the dispersion parameterisation.

These are fixed tables of the governing normative — CNCAN NSR-23, the Romanian
nuclear regulator's norm the 2021 thesis was written against — reproduced as
Tables 1-4 and 7 of the thesis. They are compile-time constants rather than data files: they
are small, they never vary between runs, and the values are read inside the
innermost loops of every field evaluation, where a keyed lookup into a data
frame costs more than the dispersion calculation itself.

The vertical dispersion tables are Hosker's, and the normative reproduces them:
[Hosker1974](@citet), IAEA-SM-181/19 — an analytic fit to the schemes of F.B.
Smith [Smith1973](@cite) and Briggs [Briggs1973](@cite). They are printed in
HPA-RPD-058 [SmithSimmonds2009](@cite) Table 3.3 and in NRPB-R91
[Clarke1979](@cite) Table 3, and both printings are asserted in the test suite.

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

"""
    DispersionScheme

Which published set of dispersion parameters a [`Site`](@ref) evaluates.

  - `DISPERSION_HOSKER` — the normative's pairing, and the default: Briggs'
    open-country `σ_y` with Hosker's `σ_z`, which carries an explicit roughness
    correction. See [`lateral_dispersion`](@ref) and
    [`vertical_dispersion`](@ref).
  - `DISPERSION_BRIGGS_OPEN_COUNTRY` — both parameters from [Briggs1973](@citet),
    as printed in the Handbook on Atmospheric Diffusion [Hanna1982](@cite)
    Table 4.5. No roughness dependence.
  - `DISPERSION_BRIGGS_URBAN` — Briggs' urban set from the same table, his fit
    to the St. Louis experiment of [McElroyPooler1968](@citet). Classes A and B
    share one curve, as do E and F.
  - `DISPERSION_EIMUTIS_KONICEK` — the analytic fit of [EimutisKonicek1972](@citet)
    to the Pasquill–Gifford curves, which the NRC codes XOQDOQ
    [Sagendorf1982](@cite) and PAVAN [Bander1982](@cite) evaluate.

Every scheme states a range of downwind distance it was fitted over;
[`validity_range`](@ref) returns it.
"""
@enum DispersionScheme::UInt8 begin
    DISPERSION_HOSKER = 1
    DISPERSION_BRIGGS_OPEN_COUNTRY = 2
    DISPERSION_BRIGGS_URBAN = 3
    DISPERSION_EIMUTIS_KONICEK = 4
end

"""
    DISPERSION_SCHEMES

The four dispersion schemes, in enumeration order.
"""
const DISPERSION_SCHEMES = (
    DISPERSION_HOSKER,
    DISPERSION_BRIGGS_OPEN_COUNTRY,
    DISPERSION_BRIGGS_URBAN,
    DISPERSION_EIMUTIS_KONICEK,
)

"""
    BriggsCoefficients

Coefficients of one of Briggs' interpolation formulae for a dispersion
parameter,

    σ = a x (1 + b x)^p

with `x` and `σ` in metres. See [`briggs_lateral_coefficients`](@ref) and
[`briggs_vertical_coefficients`](@ref).
"""
struct BriggsCoefficients
    a::Float64
    b::Float64
    p::Float64
end

# Briggs (1973), Diffusion Estimation for Small Emissions, Appendix D, as printed
# in Hanna, Briggs and Hosker (1982) Table 4.5. Briggs tabulates plume
# half-widths R = 1.25 σ; the Handbook prints the σ themselves, rounded to two
# figures, and that printing is what EPA ISC3 (Tables 1-3 and 1-4) and this
# package carry. The Handbook misprints the coefficient inside the E–F urban σ_z
# as 0.00015; Briggs' own table and ISC3 have 0.0015, and 0.0015 is used.
const _BRIGGS_OPEN_COUNTRY_LATERAL = (
    BriggsCoefficients(0.22, 1e-4, -0.5),      # A
    BriggsCoefficients(0.16, 1e-4, -0.5),      # B
    BriggsCoefficients(0.11, 1e-4, -0.5),      # C
    BriggsCoefficients(0.08, 1e-4, -0.5),      # D
    BriggsCoefficients(0.06, 1e-4, -0.5),      # E
    BriggsCoefficients(0.04, 1e-4, -0.5),      # F
)
const _BRIGGS_OPEN_COUNTRY_VERTICAL = (
    BriggsCoefficients(0.20, 0.0, 0.0),        # A
    BriggsCoefficients(0.12, 0.0, 0.0),        # B
    BriggsCoefficients(0.08, 2e-4, -0.5),      # C
    BriggsCoefficients(0.06, 1.5e-3, -0.5),    # D
    BriggsCoefficients(0.03, 3e-4, -1.0),      # E
    BriggsCoefficients(0.016, 3e-4, -1.0),     # F
)
const _BRIGGS_URBAN_LATERAL = (
    BriggsCoefficients(0.32, 4e-4, -0.5),      # A
    BriggsCoefficients(0.32, 4e-4, -0.5),      # B, as A
    BriggsCoefficients(0.22, 4e-4, -0.5),      # C
    BriggsCoefficients(0.16, 4e-4, -0.5),      # D
    BriggsCoefficients(0.11, 4e-4, -0.5),      # E
    BriggsCoefficients(0.11, 4e-4, -0.5),      # F, as E
)
const _BRIGGS_URBAN_VERTICAL = (
    BriggsCoefficients(0.24, 1e-3, 0.5),       # A
    BriggsCoefficients(0.24, 1e-3, 0.5),       # B, as A
    BriggsCoefficients(0.20, 0.0, 0.0),        # C
    BriggsCoefficients(0.14, 3e-4, -0.5),      # D
    BriggsCoefficients(0.08, 1.5e-3, -0.5),    # E
    BriggsCoefficients(0.08, 1.5e-3, -0.5),    # F, as E
)

function _briggs_tables(scheme::DispersionScheme)
    scheme == DISPERSION_BRIGGS_URBAN &&
        return (_BRIGGS_URBAN_LATERAL, _BRIGGS_URBAN_VERTICAL)
    scheme == DISPERSION_BRIGGS_OPEN_COUNTRY &&
        return (_BRIGGS_OPEN_COUNTRY_LATERAL, _BRIGGS_OPEN_COUNTRY_VERTICAL)
    throw(ArgumentError("$scheme is not one of Briggs' interpolation schemes"))
end

"""
    briggs_lateral_coefficients(class, scheme)

[`BriggsCoefficients`](@ref) of `σ_y` for the given Pasquill `class` under
`DISPERSION_BRIGGS_OPEN_COUNTRY` or `DISPERSION_BRIGGS_URBAN`.
"""
briggs_lateral_coefficients(class::PasquillClass,
    scheme::DispersionScheme,) = _briggs_tables(scheme)[1][classindex(class)]

"""
    briggs_vertical_coefficients(class, scheme)

[`BriggsCoefficients`](@ref) of `σ_z` for the given Pasquill `class` under
`DISPERSION_BRIGGS_OPEN_COUNTRY` or `DISPERSION_BRIGGS_URBAN`.
"""
briggs_vertical_coefficients(class::PasquillClass,
    scheme::DispersionScheme,) = _briggs_tables(scheme)[2][classindex(class)]

"""
    lateral_coefficient(class)

Coefficient `c₃` of the lateral dispersion parameter for the given Pasquill
class, Briggs' open-country `a`; see [`lateral_dispersion`](@ref).
"""
lateral_coefficient(class::PasquillClass) = _BRIGGS_OPEN_COUNTRY_LATERAL[classindex(class)].a

"""
    PowerLawCoefficients

Coefficients of a dispersion parameter written as

    σ = a x^b + c

with `x` and `σ` in metres, the form of the Eimutis–Konicek fit. See
[`eimutis_konicek_vertical_coefficients`](@ref).
"""
struct PowerLawCoefficients
    a::Float64
    b::Float64
    c::Float64
end

# Eimutis and Konicek (1972), Atmospheric Environment 6, 859–863, as carried in
# the DATA statements of subroutine POLYN of NRC XOQDOQ (NUREG/CR-2919) and of
# PAVAN (NUREG/CR-2858), which agree digit for digit. σ_y = a x^0.9031 at every
# distance; σ_z = a x^b + c on three ranges of distance, and the three pieces of
# each class join to within a per cent at 100 m and 1 km.
const _EIMUTIS_KONICEK_LATERAL = (0.3658, 0.2751, 0.2089, 0.1471, 0.1046, 0.0722)

"""
    EIMUTIS_KONICEK_LATERAL_EXPONENT

The exponent of distance in the Eimutis–Konicek `σ_y`, 0.9031, common to every
class.
"""
const EIMUTIS_KONICEK_LATERAL_EXPONENT = 0.9031

"""
    EIMUTIS_KONICEK_RANGES

Distances in metres at which the Eimutis–Konicek `σ_z` changes from one set of
coefficients to the next: below 100 m, from 100 m to 1 km, and beyond.
"""
const EIMUTIS_KONICEK_RANGES = (100.0, 1000.0)

const _EIMUTIS_KONICEK_VERTICAL = (
    (   # A
        PowerLawCoefficients(0.192, 0.936, 0.0),
        PowerLawCoefficients(0.00066, 1.941, 9.27),
        PowerLawCoefficients(0.00024, 2.094, -9.6),
    ),
    (   # B
        PowerLawCoefficients(0.156, 0.922, 0.0),
        PowerLawCoefficients(0.0382, 1.149, 3.3),
        PowerLawCoefficients(0.055, 1.098, 2.0),
    ),
    (   # C
        PowerLawCoefficients(0.116, 0.905, 0.0),
        PowerLawCoefficients(0.113, 0.911, 0.0),
        PowerLawCoefficients(0.113, 0.911, 0.0),
    ),
    (   # D
        PowerLawCoefficients(0.079, 0.881, 0.0),
        PowerLawCoefficients(0.222, 0.725, -1.7),
        PowerLawCoefficients(1.26, 0.516, -13.0),
    ),
    (   # E
        PowerLawCoefficients(0.063, 0.871, 0.0),
        PowerLawCoefficients(0.211, 0.678, -1.3),
        PowerLawCoefficients(6.73, 0.305, -34.0),
    ),
    (   # F
        PowerLawCoefficients(0.053, 0.814, 0.0),
        PowerLawCoefficients(0.086, 0.74, -0.35),
        PowerLawCoefficients(18.05, 0.18, -48.6),
    ),
)

"""
    eimutis_konicek_lateral_coefficient(class)

Coefficient `a` of the Eimutis–Konicek `σ_y = a x^0.9031` for the given
Pasquill `class`.
"""
eimutis_konicek_lateral_coefficient(class::PasquillClass) = _EIMUTIS_KONICEK_LATERAL[classindex(class)]

"""
    eimutis_konicek_vertical_coefficients(class, x)

[`PowerLawCoefficients`](@ref) of the Eimutis–Konicek `σ_z` for the given
Pasquill `class` on the range of distance that contains `x` metres; see
[`EIMUTIS_KONICEK_RANGES`](@ref).
"""
function eimutis_konicek_vertical_coefficients(class::PasquillClass, x::Real)
    x ≥ 0 || throw(DomainError(x, "downwind distance cannot be negative"))
    range = 1 + count(≤(x), EIMUTIS_KONICEK_RANGES)
    return _EIMUTIS_KONICEK_VERTICAL[classindex(class)][range]
end

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
profile_exponent(surface::WindProfileSurface,
    class::PasquillClass,) = _PROFILE_EXPONENT[Int(surface)][classindex(class)]

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
    TECDOC-379 [IAEA1986](@cite) §3.5.4 gives 5 × 10⁻⁸ s⁻¹ for inorganic iodine in powder snow at
    0.2 mm/h against 1.7 × 10⁻⁵ for the same species in rain, a factor of some
    340 in the same direction.
  - `WASHOUT_HTO` — the same rain columns, with snow from [Ogram1985](@citet),
    §6.0 Eq. (38):

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
equivalent, from [Ogram1985](@citet) Eq. (38).

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
