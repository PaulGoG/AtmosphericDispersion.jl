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
                   activity, precipitation = PRECIPITATION_RAIN,
                   rate = first(PRECIPITATION_RATES), washout_duration = 0,
                   release_duration = SHORT_RELEASE_REFERENCE)

Wet deposition per unit ground area in Bq/m² beneath a plume of a single
direction,

    ω_w = Λ A D exp(−y²/2Σ_y²) / (√(2π) Σ_y u)

where `Λ` is the high washout coefficient, `A` the released activity in Bq and
`D` the surviving fraction after decay and washout.

The denominator is `√(2π) Σ_y u`, which is what integrating the Gaussian plume
over the whole vertical column leaves. The 2021 code wrote `√2 π Σ_y u`; the
ratio of the two is exactly `√π`, so that expression understated wet deposition
by a factor of 1.772.
"""
function wet_deposition(
    east::Real,
    north::Real,
    site::Site,
    class::PasquillClass,
    wind_bearing::Real,
    nuclide::Nuclide;
    activity::Real,
    precipitation::PrecipitationType = PRECIPITATION_RAIN,
    rate::Real = first(PRECIPITATION_RATES),
    washout_duration::Real = 0.0,
    release_duration::Real = SHORT_RELEASE_REFERENCE,
)
    activity ≥ 0 || throw(DomainError(activity, "released activity cannot be negative"))
    x, y = plume_frame(east, north, wind_bearing)
    x > 0 || return 0.0

    Λ = washout_coefficients(precipitation, rate).high
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0
    Σy = corrected_lateral_dispersion(x, site, class; release_duration)

    surviving = decay_factor(x, u, nuclide)
    washout_duration > 0 &&
        (surviving *= wet_depletion_factor(washout_duration, precipitation, rate))

    return Λ * activity * surviving * exp(-y^2 / (2Σy^2)) / (sqrt(2π) * Σy * u)
end

"""
    wet_deposition_sector(r, site, class, nuclide;
                          activity, precipitation = PRECIPITATION_RAIN,
                          rate = first(PRECIPITATION_RATES), washout_duration = 0,
                          sectors = SectorGrid(16))

Wet deposition per unit ground area in Bq/m² at radial distance `r` metres,
averaged over a sector,

    ω_w = Λ A D / (u θ_L r)

the column activity spread over the arc the sector subtends at that distance.
"""
function wet_deposition_sector(
    r::Real,
    site::Site,
    class::PasquillClass,
    nuclide::Nuclide;
    activity::Real,
    precipitation::PrecipitationType = PRECIPITATION_RAIN,
    rate::Real = first(PRECIPITATION_RATES),
    washout_duration::Real = 0.0,
    sectors::SectorGrid = SectorGrid(16),
)
    activity ≥ 0 || throw(DomainError(activity, "released activity cannot be negative"))
    r > 0 || return 0.0

    Λ = washout_coefficients(precipitation, rate).high
    u = transport_wind_speed(site, class)
    u > 0 || return 0.0

    surviving = decay_factor(r, u, nuclide)
    washout_duration > 0 &&
        (surviving *= wet_depletion_factor(washout_duration, precipitation, rate))

    return Λ * activity * surviving / (u * sector_width(sectors) * r)
end

"""
    RESUSPENSION_COEFFICIENTS

The two amplitudes `(A, B)` in m⁻¹ of the resuspension factor, and the two
decay constants `(λ₁, λ₂)` in day⁻¹ that go with them.

The fast term describes material still loose on the surface, the slow term
material progressively fixed into it.
"""
const RESUSPENSION_COEFFICIENTS = (A = 1e-5, B = 1e-9, λ₁ = 1e-2, λ₂ = 2e-5)

"""
    ResuspensionModel

Which published resuspension correlation to use.

  - `RESUSPENSION_IAEA_SS57` — IAEA Safety Series No. 57 (1982), §3.6,
    Eq. (3.14A), the two-term form above. The default, and the one whose
    constants this package carries exactly. Safety Series 57 reached this
    package by way of reference [3] of CNCAN NSR-23; it has since been
    superseded, and every page of it carries a "no longer valid" stamp, but no
    successor restates these constants.
  - `RESUSPENSION_MAXWELL_ANSPAUGH` — Maxwell and Anspaugh, *Health Physics*
    **101** (2011), Eqs. 15/16, also adopted by NRC NUREG/CR-7270:
    `10⁻⁵e^(−0.07t) + 7×10⁻⁹e^(−0.002t) + 10⁻⁹`. It keeps the amplitudes and
    weathers seven times faster, so it falls well below Safety Series 57 from
    about a week onwards.
"""
@enum ResuspensionModel begin
    RESUSPENSION_IAEA_SS57 = 1
    RESUSPENSION_MAXWELL_ANSPAUGH = 2
end

"""
    resuspension_factor(elapsed_days)

Resuspension factor `K` in m⁻¹ at `elapsed_days` after deposition,

    K = A exp(−λ₁ t) + B exp(−λ₂ t)

It converts a surface deposition in Bq/m² into an airborne concentration in
Bq/m³. It falls by a factor of about 38 over the first year and by four orders
of magnitude over ten, as the deposited material weathers into the surface.
"""
function resuspension_factor(
    elapsed_days::Real,
    model::ResuspensionModel = RESUSPENSION_IAEA_SS57,
)
    elapsed_days ≥ 0 || throw(DomainError(elapsed_days, "elapsed time cannot be negative"))
    if model == RESUSPENSION_MAXWELL_ANSPAUGH
        return 1e-5 * exp(-0.07 * elapsed_days) + 7e-9 * exp(-0.002 * elapsed_days) + 1e-9
    end
    c = RESUSPENSION_COEFFICIENTS
    return c.A * exp(-c.λ₁ * elapsed_days) + c.B * exp(-c.λ₂ * elapsed_days)
end

"""
    resuspended_concentration(deposition, elapsed_days)

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
