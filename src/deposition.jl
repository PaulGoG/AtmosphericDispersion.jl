#=
Ground deposition and resuspension.

Deposition is expressed per unit ground area, in becquerel per square metre,
for a release of `activity` becquerel. Dry deposition is the ground-level
time-integrated concentration times a deposition velocity; wet deposition is
the washout coefficient times the whole vertical column, since rain falling
through the plume scavenges all of it, not only the part near the ground.

Both take the *high* member of their bracketing pair — the high deposition
velocity, the high washout coefficient — because more deposited material is the
conservative assumption for the ground-shine and ingestion pathways. Depletion
of the airborne plume takes the low members, for the opposite reason. See
`DepositionVelocity`.
=#

"""
    dry_deposition(concentration, nuclide)

Dry deposition per unit ground area in Bq/m², from the time-integrated
ground-level concentration in Bq·s/m³.

Takes the high deposition velocity of the nuclide.
"""
function dry_deposition(concentration::Real, nuclide::Nuclide)
    concentration ≥ 0 ||
        throw(DomainError(concentration, "concentration cannot be negative"))
    return nuclide.deposition_velocity.high * concentration
end

"""
    wet_deposition(east, north, site, class, wind_bearing, nuclide;
                   activity, washout, release_duration = SHORT_RELEASE_REFERENCE)

Wet deposition per unit ground area in Bq/m² beneath a plume of a single
direction,

    ω_w = Λ A D exp(−y²/2Σ_y²) / (√(2π) Σ_y u)

where `Λ` is the high washout coefficient of the [`WashoutEvent`](@ref)
`washout` for the nuclide's species, `A` the released activity in Bq and `D`
the surviving fraction after decay and washout.

The denominator is `√(2π) Σ_y u`, which is what integrating the Gaussian plume
over the whole vertical column leaves. The 2021 code wrote `√2 π Σ_y u`; the
ratio of the two is exactly `√π`, so that expression understated wet deposition
by a factor of 1.772.
"""
function wet_deposition(
    east::Real,
    north::Real,
    site::AbstractSite,
    class::PasquillClass,
    wind_bearing::Real,
    nuclide::Nuclide;
    activity::Real,
    washout::WashoutEvent,
    release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    activity ≥ 0 || throw(DomainError(activity, "released activity cannot be negative"))
    x, y = plume_frame(east, north, wind_bearing)
    x > 0 || return 0.0

    Λ = washout_coefficients(washout, nuclide.washout_species).high
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0
    Σy = corrected_lateral_dispersion(x, site, class; release_duration)

    surviving =
        decay_factor(x, u, nuclide) * wet_depletion_factor(washout, nuclide.washout_species)

    return Λ * activity * surviving * exp(-y^2 / (2Σy^2)) / (sqrt(2π) * Σy * u)
end

"""
    wet_deposition_sector(r, site, class, nuclide;
                          activity, washout, sectors = SectorGrid(16))

Wet deposition per unit ground area in Bq/m² at radial distance `r` metres,
averaged over a sector,

    ω_w = Λ A D / (u θ_L r)

the column activity spread over the arc the sector subtends at that distance.
"""
function wet_deposition_sector(
    r::Real,
    site::AbstractSite,
    class::PasquillClass,
    nuclide::Nuclide;
    activity::Real,
    washout::WashoutEvent,
    sectors::SectorGrid = SectorGrid(16),
)
    activity ≥ 0 || throw(DomainError(activity, "released activity cannot be negative"))
    r > 0 || return 0.0

    Λ = washout_coefficients(washout, nuclide.washout_species).high
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0

    surviving =
        decay_factor(r, u, nuclide) * wet_depletion_factor(washout, nuclide.washout_species)

    return Λ * activity * surviving / (u * sector_width(sectors) * r)
end

"""
    ResuspensionModel(; fast_amplitude, fast_rate, slow_amplitude, slow_rate, floor = 0)

A resuspension factor of the form

    K(t) = A exp(−λ₁ t) + B exp(−λ₂ t) + C

in m⁻¹, with `t` in days since deposition: amplitudes `A` and `B` in m⁻¹, rates
`λ₁` and `λ₂` in day⁻¹, and a long-term `floor` `C` in m⁻¹. The fast term
describes material still loose on the surface, the slow term material
progressively fixed into it.

See [`RESUSPENSION_IAEA_SS57`](@ref) and [`RESUSPENSION_MAXWELL_ANSPAUGH`](@ref)
for the two published parameter sets carried; any other is a call away.
"""
struct ResuspensionModel
    fast_amplitude::Float64
    fast_rate::Float64
    slow_amplitude::Float64
    slow_rate::Float64
    floor::Float64

    function ResuspensionModel(;
        fast_amplitude::Real,
        fast_rate::Real,
        slow_amplitude::Real,
        slow_rate::Real,
        floor::Real = 0.0,
    )
        all(≥(0), (fast_amplitude, fast_rate, slow_amplitude, slow_rate, floor)) || throw(
            ArgumentError("resuspension amplitudes, rates and floor cannot be negative"),
        )
        return new(fast_amplitude, fast_rate, slow_amplitude, slow_rate, floor)
    end
end

"""
    RESUSPENSION_IAEA_SS57

IAEA Safety Series No. 57 (1982), §3.6, Eq. (3.14A):
`10⁻⁵ exp(−10⁻² t) + 10⁻⁹ exp(−2×10⁻⁵ t)`. The default. Safety Series 57 reached
this package by way of reference [3] of CNCAN NSR-23; it has since been
superseded, and every page of it carries a "no longer valid" stamp, but no
successor restates these constants. It falls by a factor of about 38 over the
first year and by four orders of magnitude over ten.
"""
const RESUSPENSION_IAEA_SS57 = ResuspensionModel(;
    fast_amplitude = 1e-5,
    fast_rate = 1e-2,
    slow_amplitude = 1e-9,
    slow_rate = 2e-5,
)

"""
    RESUSPENSION_MAXWELL_ANSPAUGH

Maxwell and Anspaugh, *Health Physics* **101** (2011), Eqs. 15/16, also adopted
by NRC NUREG/CR-7270: `10⁻⁵ exp(−0.07 t) + 7×10⁻⁹ exp(−0.002 t) + 10⁻⁹`. It keeps
the fast amplitude and weathers seven times faster, so it lies below Safety
Series 57 from the first days, by a factor of 60 after a year. Its floor takes
over at about two and a half years, and from there on it is the higher of the
two, by 8 % after ten years.
"""
const RESUSPENSION_MAXWELL_ANSPAUGH = ResuspensionModel(;
    fast_amplitude = 1e-5,
    fast_rate = 0.07,
    slow_amplitude = 7e-9,
    slow_rate = 0.002,
    floor = 1e-9,
)

"""
    resuspension_factor(elapsed_days, model = RESUSPENSION_IAEA_SS57)

Resuspension factor `K` in m⁻¹ at `elapsed_days` after deposition; see
[`ResuspensionModel`](@ref). It converts a surface deposition in Bq/m² into an
airborne concentration in Bq/m³.
"""
function resuspension_factor(
    elapsed_days::Real,
    model::ResuspensionModel = RESUSPENSION_IAEA_SS57,
)
    elapsed_days ≥ 0 || throw(DomainError(elapsed_days, "elapsed time cannot be negative"))
    return model.fast_amplitude * exp(-model.fast_rate * elapsed_days) +
           model.slow_amplitude * exp(-model.slow_rate * elapsed_days) +
           model.floor
end

"""
    resuspended_concentration(deposition, elapsed_days, model = RESUSPENSION_IAEA_SS57)

Airborne concentration in Bq/m³ resuspended from a surface deposition of
`deposition` Bq/m², `elapsed_days` after it was laid down.
"""
function resuspended_concentration(
    deposition::Real,
    elapsed_days::Real,
    model::ResuspensionModel = RESUSPENSION_IAEA_SS57,
)
    deposition ≥ 0 || throw(DomainError(deposition, "deposition cannot be negative"))
    return deposition * resuspension_factor(elapsed_days, model)
end
