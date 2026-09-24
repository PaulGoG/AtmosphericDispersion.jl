# Regulatory screening tables and code outputs reproduced end to end: the
# IAEA SRS-19 tables and worked examples, and NRC XOQDOQ's own test case.
# Where the literature tests assert one constant or one formula at a time,
# these run a published calculation through the package and compare the
# printed numbers.

# IAEA Safety Reports Series No. 19 (2001): the screening tables of §3 and
# the worked examples of Annex IV. Table I tabulates the 30°-sector
# diffusion factor F of Eq. (3) for seven bands of release height, Table
# II the wake-corrected factor B of Eqs. (4)–(6) for ten bands of building
# area at H = 0, both to one significant figure at eleven distances. The
# σ_z behind them are stated in the notes to Table I: Briggs' open-country
# class D form up to 45 m, which the package evaluates, and the Jülich and
# Karlsruhe power laws E x^G above it, which enter as a prescribed σ_z.
# SRS-19's screening rule (its Fig. 8) holds F at its maximum over distance
# for every receptor nearer than that maximum, and the tables are read
# through that rule.
@testset "SRS-19 Tables I and II, and Annex IV" begin
    stack = StackSource(;
        height = 1.0,
        diameter = 1.0,
        exit_velocity = 1e-9,
        exit_density = 1.2,
        exit_temperature = 288.0,
    )
    air = Atmosphere(;
        reference_speed = 1.0,
        temperature = 288.0,
        density = 1.2,
        lapse_rate = 0.0098,
        surface = SURFACE_AGRICULTURAL,
        roughness = ROUGHNESS_PASTURE,
    )
    open_country = Site(;
        source = stack,
        atmosphere = air,
        mixing = MIXING_UNBOUNDED,        # SRS-19 has no lid
        dispersion = DISPERSION_BRIGGS_OPEN_COUNTRY,
    )
    twelve = SectorGrid(12)
    xs = (100.0, 200.0, 400.0, 800.0, 1000.0, 2000.0,
        4000.0, 8000.0, 10_000.0, 15_000.0, 20_000.0,)
    sig1(v) = round(v; sigdigits = 1)
    # The σ_z of the two elevated bands: E = 0.215, G = 0.885 for 46–80 m
    # and E = 0.265, G = 0.818 above 80 m.
    elevated_σz(x, H) = H ≤ 80 ? 0.215 * x^0.885 : 0.265 * x^0.818
    # F per unit wind speed, at one distance.
    function F(x, H)
        plume = if H ≤ 45
            PrescribedPlume(open_country; height = H, wind = 1.0)
        else
            PrescribedPlume(open_country; height = H, wind = 1.0, vertical = elevated_σz(x, H))
        end
        return dilution_extended(0.0, -x, plume, PASQUILL_D, 0.0, twelve)
    end
    # The screening rule: F held at its maximum over distance nearer than it.
    fine = 10 .^ range(2, log10(20_000); length = 400)
    function screened(H)
        profile = [F(x, H) for x in fine]
        i = argmax(profile)
        return x -> x < fine[i] ? profile[i] : F(x, H)
    end

    table_I = (
        (0.0, 5.0, (3e-3, 7e-4, 2e-4, 6e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (6.0, 15.0, (2e-3, 6e-4, 2e-4, 6e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (16.0, 25.0, (2e-4, 2e-4, 1e-4, 5e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (26.0, 35.0, (8e-5, 8e-5, 8e-5, 4e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (36.0, 45.0, (3e-5, 3e-5, 3e-5, 3e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 3e-7)),
        (46.0, 80.0, (2e-5, 2e-5, 2e-5, 2e-5, 1e-5, 4e-6, 1e-6, 3e-7, 2e-7, 1e-7, 6e-8)),
        (81.0, 200.0, (1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 5e-6, 2e-6, 5e-7, 3e-7, 1e-7, 9e-8)),
    )
    # Every cell is reproduced, to the printed figure, by a release height
    # inside its band; the table does not say which height it used.
    for (lo, hi, printed) in table_I
        tables = [screened(H) for H in range(lo, hi; length = 21)]
        for (k, x) in enumerate(xs)
            @test any(t -> sig1(t(x)) == printed[k], tables)
        end
    end

    # Table II: the same at H = 0 with Σ_z = (σ_z² + A_B/π)^(1/2), Eq. (6).
    function B(x, A)
        envelope = BuildingEnvelope([
            Building(; east = 10.0, north = 0.0, height = 10.0, frontal_area = A),
        ])
        built = Site(;
            source = stack,
            atmosphere = air,
            mixing = MIXING_UNBOUNDED,
            dispersion = DISPERSION_BRIGGS_OPEN_COUNTRY,
            buildings = envelope,
        )
        plume = PrescribedPlume(built; height = 0.0, wind = 1.0)
        return dilution_extended(0.0, -x, plume, PASQUILL_D, 0.0, twelve)
    end
    table_II = (
        (1.0, 100.0, (3e-3, 7e-4, 2e-4, 6e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (101.0, 400.0, (2e-3, 6e-4, 2e-4, 6e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (401.0, 800.0, (1e-3, 5e-4, 2e-4, 6e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (801.0, 1200.0, (9e-4, 4e-4, 2e-4, 5e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
        (1201.0, 1600.0,
            (8e-4, 3e-4, 1e-4, 5e-5, 4e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7),),
        (1601.0, 2000.0,
            (7e-4, 3e-4, 1e-4, 5e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7),),
        (2001.0, 3000.0,
            (6e-4, 3e-4, 1e-4, 5e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7),),
        (3001.0, 4000.0,
            (5e-4, 2e-4, 1e-4, 4e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7),),
        (4001.0, 6000.0,
            (4e-4, 2e-4, 9e-5, 4e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7),),
        (6001.0, Inf, (3e-4, 2e-4, 8e-5, 4e-5, 3e-5, 1e-5, 4e-6, 1e-6, 1e-6, 5e-7, 4e-7)),
    )
    # The table evaluates each band at its lower edge, the smaller area giving
    # the larger factor: there, 106 of the 110 cells round to the printed
    # figure. The four that do not sit on a rounding boundary: at 100 m for
    # 1201–1600 m² the value is 7.5 × 10⁻⁴, printed 8; and at 20 km the table
    # repeats Table I's 4 × 10⁻⁷ across every area band, where the wake term
    # is a 3 % effect and the value is 3.5 × 10⁻⁷ or below. Every cell lies
    # within the printed figure's half-unit widened by a tenth of it.
    exact = 0
    for (lo, _, printed) in table_II
        for (k, x) in enumerate(xs)
            v = B(x, lo)
            sig1(v) == printed[k] && (exact += 1)
            unit = 10.0^floor(log10(printed[k]))
            @test abs(v - printed[k]) ≤ 0.5 * unit + 0.1 * printed[k]
        end
    end
    @test exact == 106

    # Annex IV: three worked examples, ¹³¹I at 1 Bq/s, P_p = 0.25, u = 2 m/s.
    # IV-1: a 60 m stack beside 20 m buildings, a farm at 1 km. H > 2.5 H_B
    # puts it in the displacement zone; SRS-19 reads F = 10⁻⁵ m⁻² from
    # Table I and finds 1.3 × 10⁻⁶ Bq/m³.
    buildings = BuildingEnvelope([
        Building(; east = 10.0, north = 0.0, height = 20.0, frontal_area = 500.0),
    ])
    @test building_zone(60.0, 1000.0, buildings) === DISPLACEMENT_ZONE
    f = F(1000.0, 60.0)
    @test sig1(f) == 1e-5
    @test 0.25 * f / 2 ≈ 1.6e-6 rtol = 0.02       # unrounded; the example carries the table's rounding
    # IV-2: a 0.5 m vent in the side of a 500 m² building; an intake 5 m
    # along the same wall, 30/(2 · 5²) = 0.6 Bq/m³; the farm at 1 km in the
    # wake zone, B = 4 × 10⁻⁵ from Table II and 5 × 10⁻⁶ Bq/m³.
    @test building_zone(10.0, 5.0, buildings) === CAVITY_ZONE
    @test dilution_cavity_wall(5.0, 2.0; vent_diameter = 0.5) ≈ 0.6
    @test building_zone(10.0, 1000.0, buildings) === WAKE_ZONE
    b = B(1000.0, 500.0)
    @test sig1(b) == 4e-5
    @test 0.25 * b / 2 ≈ 5e-6 rtol = 0.05
    # IV-3: a 33 m stack on a 30 m building of 5000 m²; a residence 150 m
    # downwind is in the cavity, 0.25/(π · 2 · 30) = 1.3 × 10⁻³ Bq/m³, and
    # the farm at 1 km is in the wake as in IV-2.
    large = BuildingEnvelope([
        Building(; east = 10.0, north = 0.0, height = 30.0, frontal_area = 5000.0),
    ])
    @test building_zone(33.0, 150.0, large) === CAVITY_ZONE
    @test building_zone(33.0, 1000.0, large) === WAKE_ZONE
    @test 0.25 * dilution_cavity(2.0, 30.0) ≈ 1.3e-3 rtol = 0.02
end

# NRC XOQDOQ, NUREG/CR-2919 (1982), Test Case 2: a continuous elevated release
# evaluated by the code itself, Appendix B (input) and Appendix C (output). A
# 10 m/s jet from a 2 m stack at 45 m with no heat; a joint frequency
# distribution of 100 hours in three directions, five wind-speed classes and
# stability classes C to G; terrain rising to 200 m at 10 km; no recirculation
# correction. The undepleted annual-average χ/Q is printed at 22 distances to
# four figures. Everything XOQDOQ does on the way is stated in the report and
# reproduced here from it: the wind adjusted from 10 m to the release height
# by (45/10)^0.25 in classes A–D and ^0.5 in E–G (§4.17); Briggs' momentum
# rise, Eqs. 16–20, with the downwash of Eq. 17 and the exponents as the
# FORTRAN rounds them; the terrain interpolated linearly and the effective
# height floored at zero (Eq. 15, subroutine HEIGHT); the Eimutis–Konicek σ
# capped at 1000 m and class G taken as σ_F²/σ_E (subroutine POLYN); and the
# sector-averaged kernel with 2.032, which is the package's own. The plume
# enters the package prescribed: its height, wind and σ_z per class.
@testset "XOQDOQ Test Case 2" begin
    mile = 1609.347219                    # the code's own constant
    miles = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
        5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,)
    # Sector S, the wind from N at 25 % of the hours.
    printed = (2.372e-7, 6.713e-7, 1.045e-6, 1.377e-6, 1.670e-6, 1.322e-6, 9.637e-7,
        7.438e-7, 5.969e-7, 4.929e-7, 4.166e-7, 3.588e-7, 2.037e-7, 1.374e-7,
        7.982e-8, 5.484e-8, 4.111e-8, 3.254e-8, 2.673e-8, 2.257e-8, 1.945e-8,
        1.703e-8,)
    stack_height, exit_velocity, diameter = 45.0, 10.0, 2.0
    wind_height, release_height_wind = 10.0, 45.0
    class_speeds = (0.5, 1.5, 3.0, 6.0, 12.0)         # class midpoints at 10 m
    terrain = ((100.0, 0.0), (800.0, 16.0), (10_000.0, 200.0))   # (distance, height)
    # Hours in each cell of the joint frequency distribution for this sector,
    # classes C to G by wind class, out of 100.
    hours = 1.0
    total_hours = 100.0

    # Subroutine HEIGHT: linear interpolation between the terrain points, from
    # plant grade at the site, extrapolated along the last segment beyond.
    function terrain_height(x)
        m = something(findfirst(p -> first(p) ≥ x, terrain), length(terrain))
        d₀, h₀ = m == 1 ? (0.0, 0.0) : terrain[m - 1]
        d₁, h₁ = terrain[m]
        return h₀ + (h₁ - h₀) * (x - d₀) / (d₁ - d₀)
    end
    # Subroutine RISE with no heat: Briggs' momentum rise, the stable limits in
    # E to G with the code's own stability parameters, the neutral limit
    # 3 w₀D/u, and Gifford's downwash 3(1.5 − w₀/u)D below w₀ = 1.5u.
    stability = Dict(5 => 0.000875, 6 => 0.00175, 7 => 0.00245)
    function rise(x, u, class)
        Δh = 1.44 * (exit_velocity / u)^0.667 * (x / diameter)^0.333 * diameter
        if class > 4
            Fₘ = exit_velocity^2 * diameter^2 / 4
            S = stability[class]
            Δh = min(Δh, 4 * (Fₘ / S)^0.25, 1.5 * (Fₘ / u)^0.333 / S^0.1667)
        end
        Δh = min(Δh, 3 * exit_velocity * diameter / u)
        exit_velocity < 1.5u && (Δh -= 3 * (1.5 - exit_velocity / u) * diameter)
        return Δh
    end
    # Subroutine POLYN: Eimutis–Konicek, capped at 1000 m, G from E and F.
    function σz(x, class)
        class ≤ 6 &&
            return min(eimutis_konicek_vertical_dispersion(x, PASQUILL_CLASSES[class]), 1000.0)
        e = eimutis_konicek_vertical_dispersion(x, PASQUILL_E)
        f = eimutis_konicek_vertical_dispersion(x, PASQUILL_F)
        return min(f^2 / e, 1000.0)
    end

    stack = StackSource(;
        height = stack_height,
        diameter,
        exit_velocity,
        exit_density = 1.2,
        exit_temperature = 288.0,
    )
    air = Atmosphere(;
        reference_speed = 1.0,
        temperature = 288.0,
        density = 1.2,
        lapse_rate = 0.0098,
        surface = SURFACE_AGRICULTURAL,
        roughness = ROUGHNESS_PASTURE,
    )
    site = Site(;
        source = stack,
        atmosphere = air,
        mixing = MIXING_UNBOUNDED,          # XOQDOQ has no lid
        dispersion = DISPERSION_EIMUTIS_KONICEK,
    )
    sixteen = SectorGrid(16)
    # Subroutine ANNUAL: the frequency-weighted sum of the sector-averaged
    # kernel over the wind and stability classes, the effective height per
    # cell from Eq. 15 and the plume dropped where h_e/σ_z ≥ 15.
    function annual(x)
        total = 0.0
        for class in 3:7, i in 1:5

            exponent = class ≤ 4 ? 0.25 : 0.5
            u = class_speeds[i] * (release_height_wind / wind_height)^exponent
            h = max(stack_height + rise(x, u, class) - terrain_height(x), 0.0)
            σ = σz(x, class)
            h / σ ≥ 15 && continue
            label = class ≤ 6 ? PASQUILL_CLASSES[class] : PASQUILL_F   # no lid: the class label is inert
            plume = PrescribedPlume(site; height = h, wind = u, vertical = σ)
            total += hours / total_hours *
                     dilution_extended(0.0, -x, plume, label, 0.0, sixteen)
        end
        return total
    end

    # All 22 distances to within 0.05 %, the precision of the printed figures.
    for (k, r) in enumerate(miles)
        @test annual(r * mile) ≈ printed[k] rtol = 5e-4
    end
    # The sector receiving 50 % of the hours prints exactly twice these values,
    # and the third sector the same: 4.744 × 10⁻⁷ and 2.372 × 10⁻⁷ at 0.25 mile,
    # 3.406 × 10⁻⁸ and 1.703 × 10⁻⁸ at 50 miles.
    @test 2 * printed[1] ≈ 4.744e-7 rtol = 1e-3
    @test 2 * printed[end] ≈ 3.406e-8 rtol = 1e-3
    # The plume is at ground level from the terrain onwards: at 5 miles the
    # ridge stands above the highest effective height, that of the lightest
    # wind in class C.
    lightest = class_speeds[1] * (release_height_wind / wind_height)^0.25
    @test terrain_height(5 * mile) > stack_height + rise(5 * mile, lightest, 3)
end
