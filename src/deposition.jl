#=
Ground deposition and resuspension.

Deposition is expressed per unit ground area, in becquerel per square metre,
for a release of `activity` becquerel. Dry deposition is the ground-level
time-integrated concentration times a deposition velocity; wet deposition is
the washout coefficient times the whole vertical column, since rain falling
through the plume scavenges all of it.

Both take the high member of their bracketing pair (the high deposition
velocity and the high washout coefficient) because more deposited material is
the conservative assumption for the ground-shine and ingestion pathways.
Depletion of the airborne plume takes the low members, for the opposite reason.
See `DepositionVelocity`. Wet deposition depletes its own plume with the same
high coefficient, so that the activity it deposits cannot exceed the activity
released during the precipitation.
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

    ω_w = Λ f A D exp(−Λ x/u) exp(−y²/2Σ_y²) / (√(2π) Σ_y u)

where `Λ` is the high washout coefficient of the [`WashoutEvent`](@ref)
`washout` for the nuclide's species, `f` its [`precipitation_fraction`](@ref),
`A` the released activity in Bq and `D` the surviving fraction after decay. The
column is depleted by the same `Λ` over the travel time `x/u`. Beyond the
distance `u t_w` that material travels while the precipitation of duration
`t_w` lasts, it is no longer under precipitation and the deposit is zero. The
deposit integrated over the ground out to `x` is then
`f A (1 − exp(−Λ min(x/u, t_w)))` (with no decay), never more than the activity
released during the precipitation, and nothing for an event of zero duration.

CNCAN NSR-23 [CNCAN2004](@cite) Eq. (38) multiplies the high coefficient by the
whole release depleted with the low coefficient over the duration of the event,
which the 2021 code followed; that deposit does not vanish without
precipitation and its integral over the ground can exceed the release.

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
    x / u ≤ washout.duration || return 0.0
    Σy = corrected_lateral_dispersion(x, site, class; release_duration)

    column = precipitation_fraction(washout) * activity * decay_factor(x, u, nuclide) *
             exp(-Λ * x / u)
    return Λ * column * exp(-y^2 / (2Σy^2)) / (sqrt(2π) * Σy * u)
end

"""
    wet_deposition_sector(r, site, class, nuclide;
                          activity, washout, sectors = SectorGrid(16))

Wet deposition per unit ground area in Bq/m² at radial distance `r` metres,
averaged over a sector,

    ω_w = Λ f A D exp(−Λ r/u) / (u θ_L r)

the column activity spread over the arc the sector subtends at that distance,
with `f`, `D`, the depletion of the column and the zero beyond `u t_w` as in
[`wet_deposition`](@ref). NSR-23 Eq. (39) is the same form with the whole
release, depleted over the duration of the event.
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
    r / u ≤ washout.duration || return 0.0

    column = precipitation_fraction(washout) * activity * decay_factor(r, u, nuclide) *
             exp(-Λ * r / u)
    return Λ * column / (u * sector_width(sectors) * r)
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
for the two published parameter sets carried.
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

IAEA Safety Series No. 57 [IAEA1982](@cite), §3.6, Eq. (3.14A):
`10⁻⁵ exp(−10⁻² t) + 10⁻⁹ exp(−2×10⁻⁵ t)`. It is the default. Safety Series 57
was adopted here as reference [3] of CNCAN NSR-23; it has since been
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

[MaxwellAnspaugh2011](@citet), Eqs. 15/16, also adopted by NRC NUREG/CR-7270
[Bixler2022](@cite): `10⁻⁵ exp(−0.07 t) + 7×10⁻⁹ exp(−0.002 t) + 10⁻⁹`. It keeps
the fast amplitude and weathers seven times faster, so it lies below Safety
Series 57 from the first days, by a factor of 60 after a year. From about two
and a half years, where its floor becomes the dominant term, it is the higher of
the two, by 8 % after ten years.
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
