# Validation against the published literature, as distinct from the analytic
# self-consistency above. Several of the parameterisations this code inherits
# from the normative turn out to be standard published schemes, and where
# they are, the identity is asserted rather than described.
@testset "Literature validation" begin
    # The vertical dispersion scheme is Hosker's fit to F.B. Smith (1972)
    # and Briggs (1973), published as IAEA-SM-181/19 (1974) and printed in
    # Smith and Simmonds (eds.), HPA-RPD-058, Health Protection Agency
    # (2009), Table 3.3, and in Clarke, NRPB-R91 (1979), Table 3. Both
    # printings agree digit for digit, and the ORNL codes that implement it
    # (ORNL-5913 Table 6, ORNL/TM-6874 p. 25) carry the same numbers.
    #
    # Transcribing a table is exactly where a package of tabulated constants
    # fails, so the whole table is asserted rather than sampled.
    @testset "vertical dispersion coefficients are Hosker's" begin
        # σ_z = a x^b / (1 + c x^d), x and σ_z in metres.
        shape = (
            (0.112, 1.06, 5.38e-4, 0.815),    # A
            (0.130, 0.950, 6.52e-4, 0.750),   # B
            (0.112, 0.920, 9.05e-4, 0.718),   # C
            (0.098, 0.889, 1.35e-3, 0.688),   # D
            (0.0609, 0.895, 1.96e-3, 0.684),  # E
            (0.0638, 0.783, 1.36e-3, 0.672),  # F
        )
        for (i, class) in enumerate(PASQUILL_CLASSES)
            c = vertical_shape_coefficients(class)
            @test (c.a₁, c.b₁, c.a₂, c.b₂) == shape[i]
        end

        # F(z₀,x) = ln(f x^g [1 + {h x^j}⁻¹]) for z₀ > 0.1 m,
        # F(z₀,x) = ln(f x^g [1 + h x^j]⁻¹)   for z₀ ≤ 0.1 m.
        correction = (
            (0.01, 1.56, 0.0480, 6.25e-4, 0.45),
            (0.04, 2.02, 0.0269, 7.76e-4, 0.37),
            (0.10, 2.72, 0.0, 0.0, 0.0),
            (0.40, 5.16, -0.098, 18.6, -0.225),
            (1.00, 7.37, -0.0957, 4.29e3, -0.60),
            (4.00, 11.7, -0.128, 4.59e4, -0.78),
        )
        for (i, roughness) in enumerate(ROUGHNESS_CLASSES)
            c = roughness_coefficients(roughness)
            @test (c.z₀, c.c₁, c.d₁, c.c₂, c.d₂) == correction[i]
        end

        # z₀ = 10 cm is the reference the fit is normalised on, so its
        # correction is ln(2.72) — unity to within the rounding of e.
        @test roughness_correction(1000.0, ROUGHNESS_PASTURE) ≈ 1 atol = 7e-4
    end

    # The sector-averaged long-term form is a regulatory equation, stated
    # identically by NRC Regulatory Guide 1.111 Rev. 1 (1977) Eq. (3) and
    # XOQDOQ, NUREG/CR-2919 (1982) Eq. (1), by IAEA Safety Reports Series
    # No. 19 (2001) Eq. (V-2), and by the German AVV zu §47 StrlSchV (2012)
    # Eq. (4.4). RG 1.111 writes the constant as 2.032 and states in words
    # that it is √(2/π) divided by a 22.5° sector in radians; SRS-19 works
    # in twelve sectors, where the same constant is 1.5238.
    @testset "the sector constant is the published regulatory value" begin
        # SRS-19 Eq. (3) prints the twelve-sector constant as 12/√(2π³),
        # which is a closed form rather than a rounded decimal; the
        # sixteen-sector analogue is RG 1.111's 2.032.
        @test sqrt(2 / π) * 16 / (2π) ≈ 16 / sqrt(2π^3) rtol = 1e-14
        @test sqrt(2 / π) * 12 / (2π) ≈ 12 / sqrt(2π^3) rtol = 1e-14
        @test sqrt(2 / π) * 16 / (2π) ≈ 2.032 rtol = 2e-4
        @test sqrt(2 / π) * 12 / (2π) ≈ 1.5238 rtol = 2e-4

        stack = StackSource(;
            height = 50.3,
            diameter = 2.33,
            exit_velocity = 10.0,
            exit_density = 0.6,
            exit_temperature = 324.0,
        )
        air = Atmosphere(;
            reference_speed = 4.0,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = 2e-2,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        # RG 1.111 Eq. (3) and SRS-19 Eq. (V-2) are written without a lid.
        site = Site(; source = stack, atmosphere = air, mixing = MIXING_UNBOUNDED)
        distances = (500.0, 2000.0, 10_000.0)

        for (n, published) in ((16, 2.032), (12, 1.5238))
            sectors = SectorGrid(n)
            for class in PASQUILL_CLASSES, r in distances
                H = effective_height(r, site, class)
                Σz = corrected_vertical_dispersion(r, site, class)
                u = transport_wind_speed(site, class)
                regulatory = published * exp(-H^2 / (2Σz^2)) / (Σz * u * r)
                @test dilution_extended(0.0, -r, site, class, 0.0, sectors) ≈ regulatory rtol =
                    2e-4
            end
        end
    end

    # IAEA Safety Series No. 57, Generic Models and Parameters for Assessing
    # the Environmental Transfer of Radionuclides from Routine Releases,
    # STI/PUB/611, Vienna (1982), §3.6, Eq. (3.14A). Read in the primary
    # document, which states verbatim: "Reference values for the terms in
    # Eq. (3.14A) are 10⁻⁵ m⁻¹, 10⁻⁹ m⁻¹ for A and B respectively and
    # 1 × 10⁻² d⁻¹ and 2 × 10⁻⁵ d⁻¹ for λ₁ and λ₂ respectively."
    #
    # Safety Series 57 has been superseded — every page of it now carries a
    # "no longer valid" stamp — but it is the source these constants came
    # from, by way of reference [3] of CNCAN NSR-23, and no successor
    # restates them. See RESUSPENSION_MAXWELL_ANSPAUGH for the modern form.
    @testset "resuspension is IAEA Safety Series 57" begin
        c = RESUSPENSION_COEFFICIENTS
        @test c.A == 1e-5
        @test c.B == 1e-9
        @test c.λ₁ == 1e-2
        @test c.λ₂ == 2e-5

        # The same paragraph brackets them across the literature.
        @test 1e-6 ≤ c.A ≤ 1e-4
        @test 1e-10 ≤ c.B ≤ 1e-8
        # "of the order of weeks" and "in the range 50 to 100 years"
        @test 14 ≤ log(2) / c.λ₁ ≤ 120
        @test 50 ≤ log(2) / c.λ₂ / 365.25 ≤ 100

        @test resuspension_factor(0) ≈ c.A + c.B
        @test resuspension_factor(1e9) ≈ 0 atol = 1e-12
    end

    # The washout table follows the published intensity dependence and
    # brackets the published amplitudes; it does not reproduce a single
    # tabulation, and is asserted as what it is.
    @testset "washout follows the published intensity law" begin
        # Λ ∝ J^0.75, from Slinn (1977) via NRPB-R322 (ADMLC, 2001) §3.1.1:
        # "a net dependence of scavenging coefficient on J^0.75".
        for precipitation in (PRECIPITATION_RAIN, PRECIPITATION_SNOW),
            field in (:low, :high)

            Λ = [
                getfield(washout_coefficients(precipitation, j), field) for
                j in PRECIPITATION_RATES
            ]
            x = log.(collect(PRECIPITATION_RATES))
            y = log.(Λ)
            x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
            p = sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
            @test 0.65 ≤ p ≤ 0.80
        end

        # At 1 mm/h the German AVV zu §47 StrlSchV (2012), Anhang 7
        # Tabelle 3, gives 7e-5 s⁻¹ for aerosols and elemental iodine and
        # 3.5e-5 for tritiated water. Both lie inside this table's range.
        w = washout_coefficients(PRECIPITATION_RAIN, 1.0)
        @test w.low ≤ 7.0e-5 ≤ w.high
        @test w.low ≤ 3.5e-5 ≤ w.high

        # IAEA Safety Series 57 Table II gives the washout coefficient as
        # Λ = a·I with a = 1.6e-4 h(mm·s)⁻¹ for particulates and 1.1e-4 for
        # elemental iodine, so 1.6e-4 and 1.1e-4 s⁻¹ at 1 mm/h. Both fall
        # inside both species rows of the normative's rain bracket.
        for species in (WASHOUT_TRITIUM_IODINE, WASHOUT_OTHER_NUCLIDES)
            r = washout_coefficients(PRECIPITATION_RAIN, 1.0, WASHOUT_NORMATIVE, species)
            @test r.low ≤ 1.6e-4 ≤ r.high
            @test r.low ≤ 1.1e-4 ≤ r.high
        end

        # Safety Series 57 makes the dependence linear in intensity, where
        # the normative's table fits 0.75 and NRPB-R322 publishes 0.75.
        # NRPB-R157 §D3.4 brackets the exponent at 0.5 to 1.0, which
        # contains both positions, so neither is asserted against the other.
        @test 0.5 ≤ 0.75 ≤ 1.0
        @test 0.5 ≤ 1.0 ≤ 1.0

        # The snow columns are the rain columns scaled down by a constant —
        # a particle-scavenging suppression, and not a measurement.
        for (field, factor) in ((:low, 100), (:high, 500))
            ratios = [
                getfield(washout_coefficients(PRECIPITATION_RAIN, j), field) /
                getfield(washout_coefficients(PRECIPITATION_SNOW, j), field) for
                j in PRECIPITATION_RATES
            ]
            @test all(r -> isapprox(r, factor; rtol = 0.2), ratios)
        end
    end

    # Turner, Workbook of Atmospheric Dispersion Estimates, PHS 999-AP-26
    # (rev. 1970), Problem 9 and Table 7-4: the concentration profile with
    # height from the ground to 450 m, at x = 1 km on the plume axis.
    #
    # Q = 151 g/s, H = 150 m, u = 4 m/s, σ_y = 157 m, σ_z = 110 m. This is
    # the only published table found that exercises *both* reflection terms
    # at arbitrary receptor height; everything else evaluates at z = 0,
    # where the two coincide and a sign error in either would cancel.
    @testset "vertical profile against Turner Table 7-4" begin
        Q, H, u, σy, σz = 151.0, 150.0, 4.0, 157.0, 110.0
        # z in m, published χ in g/m³
        published = (
            (0, 2.78e-4),
            (30, 2.85e-4),
            (60, 3.06e-4),
            (90, 3.34e-4),
            (120, 3.55e-4),
            (150, 3.58e-4),
            (180, 3.41e-4),
            (210, 3.03e-4),
            (240, 2.51e-4),
            (270, 1.94e-4),
            (300, 1.39e-4),
            (330, 9.14e-5),
            (360, 5.64e-5),
            (390, 3.26e-5),
            (420, 1.75e-5),
            (450, 8.40e-6),
        )

        stack = StackSource(;
            height = H,
            diameter = 1.0,
            exit_velocity = 1e-9,
            exit_density = 1.2,
            exit_temperature = 288.0,
        )
        air = Atmosphere(;
            reference_speed = 4.0,
            temperature = 288.0,
            density = 1.2,
            lapse_rate = 0.0098,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        site = Site(;
            source = stack,
            atmosphere = air,
            fixed_height = H,
            fixed_wind = u,
            fixed_lateral = σy,
            fixed_vertical = σz,
            mixing = MIXING_UNBOUNDED,        # Turner has no lid
        )

        worst = 0.0
        for (z, χ) in published
            computed =
                Q * dilution_instantaneous(0.0, -1000.0, Float64(z), site, PASQUILL_D, 0.0)
            @test computed ≈ χ rtol = 0.025
            worst = max(worst, abs(computed / χ - 1))
        end
        # The residual is Turner's own rounding: his intermediate columns
        # carry three significant figures and the exponentials are summed
        # from those.
        @test worst < 0.025

        # The profile peaks at plume height and is symmetric about it only
        # in the first term; the ground reflection is what lifts z = 0 above
        # the pure Gaussian and what makes the peak sit slightly below H.
        χs = [
            Q * dilution_instantaneous(0.0, -1000.0, Float64(z), site, PASQUILL_D, 0.0)
            for (z, _) in published
        ]
        @test argmax(χs) == 6                       # z = 150 m, the release height
        @test χs[1] > Q * exp(-0.5 * (H / σz)^2) / (2π * σy * σz * u)

        # Turner prints the prefactor as 3.5e-5 g/m³ where his own table
        # requires 3.5e-4: 151/(2π·157·110·4) = 3.479e-4. A second typo in a
        # published source, after HPA Table 3.7's non-monotonic cell.
        @test Q / (2π * σy * σz * u) ≈ 3.479e-4 rtol = 1e-3
        @test published[1][2] / 0.794 ≈ 3.5e-4 rtol = 0.01
    end

    # HPA-RPD-058 Table 3.7, "Fractions of material remaining in the plume
    # due to dry deposition for a deposition velocity of 10⁻² m s⁻¹":
    # 3 effective release heights × 7 stability categories × 8 distances,
    # with the wind speed at stack height printed for each row.
    #
    # This is the only end-to-end benchmark in the suite. It exercises σ_z,
    # the mixing layer and the depletion integral together against published
    # numbers, where everything else here checks one piece at a time.
    # Category G is omitted: this package has no category G.
    @testset "depletion against HPA-RPD-058 Table 3.7" begin
        rows = (
            (30.0, PASQUILL_A, 1.32, (0.96, 0.94, 0.92, 0.89, 0.86, 0.81, 0.68, 0.51)),
            (30.0, PASQUILL_B, 2.65, (0.98, 0.96, 0.94, 0.90, 0.87, 0.83, 0.73, 0.59)),
            (30.0, PASQUILL_C, 6.62, (0.99, 0.98, 0.97, 0.95, 0.93, 0.90, 0.85, 0.78)),
            (30.0, PASQUILL_D, 6.62, (1.00, 0.98, 0.97, 0.93, 0.90, 0.86, 0.78, 0.69)),
            (30.0, PASQUILL_E, 3.97, (1.00, 0.98, 0.94, 0.86, 0.78, 0.68, 0.51, 0.36)),
            (30.0, PASQUILL_F, 2.65, (1.00, 0.99, 0.95, 0.79, 0.61, 0.40, 0.13, 0.19)),
            (70.0, PASQUILL_A, 1.64, (0.99, 0.97, 0.95, 0.93, 0.90, 0.86, 0.75, 0.59)),
            (70.0, PASQUILL_B, 3.28, (1.00, 0.99, 0.97, 0.94, 0.92, 0.88, 0.80, 0.67)),
            (70.0, PASQUILL_C, 8.21, (1.00, 1.00, 0.99, 0.97, 0.96, 0.93, 0.89, 0.83)),
            (70.0, PASQUILL_D, 8.21, (1.00, 1.00, 0.99, 0.97, 0.94, 0.91, 0.84, 0.77)),
            (70.0, PASQUILL_E, 4.93, (1.00, 1.00, 0.99, 0.95, 0.89, 0.81, 0.65, 0.49)),
            (70.0, PASQUILL_F, 3.28, (1.00, 1.00, 1.00, 0.98, 0.89, 0.68, 0.28, 0.060)),
            (100.0, PASQUILL_A, 1.80, (0.99, 0.98, 0.97, 0.94, 0.92, 0.88, 0.77, 0.62)),
            (100.0, PASQUILL_B, 3.60, (1.00, 0.99, 0.98, 0.96, 0.93, 0.90, 0.82, 0.70)),
            (100.0, PASQUILL_C, 8.99, (1.00, 1.00, 0.99, 0.98, 0.97, 0.95, 0.91, 0.85)),
            (100.0, PASQUILL_D, 8.99, (1.00, 1.00, 1.00, 0.98, 0.96, 0.93, 0.87, 0.80)),
            (100.0, PASQUILL_E, 5.40, (1.00, 1.00, 1.00, 0.98, 0.94, 0.85, 0.72, 0.55)),
            (100.0, PASQUILL_F, 3.60, (1.00, 1.00, 1.00, 1.00, 0.94, 0.76, 0.33, 0.083)),
        )
        xs = (500.0, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5)

        marker = Nuclide(;
            name = "table 3.7 marker",
            decay_constant = 0.0,
            deposition_velocity = DepositionVelocity(1e-2, 1e-2),
        )
        air = Atmosphere(;
            reference_speed = 4.0,
            temperature = 288.0,
            density = 1.2,
            lapse_rate = 0.0098,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        stack = StackSource(;
            height = 30.0,
            diameter = 1.0,
            exit_velocity = 1e-9,
            exit_density = 1.2,
            exit_temperature = 288.0,
        )

        # The one cell that is wrong in the source: category F at 30 m runs
        # 0.13 at 50 km and 0.19 at 100 km. A depletion factor cannot rise
        # with distance — material already deposited does not return — so
        # the row is non-monotonic and cannot be right. This package gives
        # 0.0188 there, which is both monotone and a decimal point away from
        # the printed 0.19.
        @test rows[6][4][8] > rows[6][4][7]        # the table, not monotone

        worst = 0.0
        for (H, class, u, published) in rows
            site =
                Site(; source = stack, atmosphere = air, fixed_height = H, fixed_wind = u)
            for (k, x) in enumerate(xs)
                (H == 30.0 && class == PASQUILL_F && k == 8) && continue
                f = dry_depletion_factor(x, site, class, marker)
                @test f ≈ published[k] atol = 0.03
                worst = max(worst, abs(f - published[k]))
            end
            # every row is monotone here, including the one that is not in
            # the table
            factors = [dry_depletion_factor(x, site, class, marker) for x in xs]
            @test issorted(factors; rev = true)
        end
        @test worst < 0.03

        # The typo cell, computed rather than read
        typo_site =
            Site(; source = stack, atmosphere = air, fixed_height = 30.0, fixed_wind = 2.65)
        @test dry_depletion_factor(1e5, typo_site, PASQUILL_F, marker) ≈ 0.019 atol = 0.002
    end

    # The notes to CNCAN NSR-23 Table 6 give the deposition velocity of
    # HTO as 0.4-0.8e-2 m/s and of HT as an order of magnitude lower,
    # 0.04-0.05e-2, attributed to Murphy, Health Physics 65(6), 1993.
    @testset "tritium deposition velocities are NSR-23 Table 6" begin
        hto = TRITIATED_WATER.deposition_velocity
        ht = TRITIUM_GAS.deposition_velocity
        @test hto.low == 0.4e-2
        @test hto.high == 0.8e-2
        @test ht.low == 0.04e-2
        @test ht.high == 0.05e-2

        # "an order of magnitude lower" — the norm's own words
        @test hto.low / ht.low == 10
        @test 10 ≤ hto.high / ht.high ≤ 20

        # The lower HTO bound is close to what other sources put it at:
        # the MACCS2 default 0.5 cm/s, the Savannah River measurement 0.42,
        # the AECL range 0.392-0.444.
        for v in (0.5e-2, 0.42e-2, 0.444e-2)
            @test hto.low ≤ v ≤ hto.high
        end
        # AECL's lower end sits just outside, 2 % below the norm's 0.4 —
        # worth pinning rather than rounding away, since it says the norm's
        # lower bound is at the edge of the measurements and not inside them.
        @test 0.392e-2 < hto.low
        @test hto.low / 0.392e-2 ≈ 1.02 rtol = 0.01
    end


    # The building-wake coefficient. IAEA Safety Reports Series No. 19
    # (2001) Eq. (6) writes the wake-broadened vertical parameter as
    # (σ_z² + A_B/π)^(1/2), and the German AVV zu §47 StrlSchV (2012)
    # Eqs. (4.31)/(4.32) as √(σ² + I_G²/π): a coefficient of one in both.
    @testset "wake coefficient is the published one" begin
        @test DEFAULT_WAKE_COEFFICIENT == 1.0
        @test NORMATIVE_WAKE_COEFFICIENT == 1.5

        b = Building(; east = 25.0, north = 0.0, height = 60.0, frontal_area = 3600.0)
        iaea = BuildingEnvelope([b])
        normative = BuildingEnvelope([b]; wake_coefficient = NORMATIVE_WAKE_COEFFICIENT)
        A = equivalent_area(iaea)
        for σ in (10.0, 50.0, 200.0)
            # released inside the cavity, where the correction is undiluted
            @test wake_broadened(σ, 0.0, iaea) ≈ sqrt(σ^2 + A / π)
            @test wake_broadened(σ, 0.0, normative) ≈ sqrt(σ^2 + 1.5A / π)
            @test wake_broadened(σ, 0.0, iaea) < wake_broadened(σ, 0.0, normative)
        end
    end

    # Ogram, Precipitation Scavenging of Tritiated Water Vapour (HTO),
    # Ontario Hydro Research Division 85-233-K (1985), §6.0 Eq. (38), with
    # its Table V at 0.5, 1 and 2 mm/h. Snow scavenging of HTO is isotopic
    # exchange at the crystal surface, so the particle suppression the
    # normative applies is the wrong physics for it.
    @testset "HTO snow washout is Ogram 1985" begin
        @test ogram_snow_washout(0.5) ≈ 2.9e-4 rtol = 0.02
        @test ogram_snow_washout(1.0) ≈ 4.2e-4 rtol = 0.02
        @test ogram_snow_washout(2.0) ≈ 6.1e-4 rtol = 0.02
        @test ogram_snow_washout(1.0) ≈ 1.2e-4 + 3.0e-4

        # Ogram's measurement lands on NSR-23's own other-nuclides snow
        # lower limit, a row it does not cite Ogram for — which is the
        # evidence that the tritium row's snow column is the anomaly.
        for r in PRECIPITATION_RATES
            other = washout_coefficients(
                PRECIPITATION_SNOW,
                r,
                WASHOUT_NORMATIVE,
                WASHOUT_OTHER_NUCLIDES,
            )
            @test 0.8 < ogram_snow_washout(r) / other.low < 1.2
        end

        # NSR-23 Table 7 has two species rows; the other-nuclides snow
        # values run three to four orders of magnitude above the tritium
        # row's, in the opposite direction to the particle suppression.
        for r in PRECIPITATION_RATES
            trit = washout_coefficients(
                PRECIPITATION_SNOW,
                r,
                WASHOUT_NORMATIVE,
                WASHOUT_TRITIUM_IODINE,
            )
            other = washout_coefficients(
                PRECIPITATION_SNOW,
                r,
                WASHOUT_NORMATIVE,
                WASHOUT_OTHER_NUCLIDES,
            )
            @test other.low / trit.low > 1e3
        end

        # rain is untouched by the choice; snow differs by about 10³
        for r in PRECIPITATION_RATES
            @test washout_coefficients(PRECIPITATION_RAIN, r, WASHOUT_HTO) ==
                  washout_coefficients(PRECIPITATION_RAIN, r, WASHOUT_NORMATIVE)
            ratio =
                washout_coefficients(PRECIPITATION_SNOW, r, WASHOUT_HTO).high /
                washout_coefficients(PRECIPITATION_SNOW, r, WASHOUT_NORMATIVE).high
            @test 900 < ratio < 1500
        end
    end

    # The combined momentum-and-buoyancy law must reduce to the pure
    # buoyancy law when the momentum flux vanishes. Briggs writes its
    # buoyancy term as 3F x²/(2β²u³) with β = 0.6, a denominator of
    # 2β² = 0.72, and that is now the default; the 2021 code had 0.5, which
    # overshoots the two-thirds law by 13.6 % and is kept as THESIS_RISE.
    @testset "combined rise reduces to the two-thirds law" begin
        # With c = 2β² = 0.72 the combined law implies a two-thirds
        # coefficient of (3/0.72)^(1/3) = 1.60915, against the 1.6 the
        # literature rounds it to and that `buoyant_rise` uses. The two
        # therefore agree to 0.6 %, not exactly — the same rounding the
        # neutral final rise shows as 21.425 against a published 21.4.
        @test (3 / BRIGGS_RISE.combined_buoyancy)^(1 / 3) ≈ 1.60915 rtol = 1e-4

        checked = 0
        for F in (5.0, 50.0, 500.0), u in (2.0, 5.0, 9.0), x in (50.0, 100.0, 300.0)
            two_thirds = 1.6 * F^(1 / 3) * x^(2 / 3) / u
            # The combined law is capped at the sum of the two final rises,
            # and the cap differs between the two coefficient sets, so only
            # assert where neither is binding.
            cap(r) =
                final_momentum_rise(1e-14, 10.0, 2.0, u, -1e-6, r) +
                final_buoyant_rise(F, u, -1e-6, r)
            briggs = combined_rise(x, F, 1e-14, 10.0, u, -1e-6, 2.0, BRIGGS_RISE)
            thesis = combined_rise(x, F, 1e-14, 10.0, u, -1e-6, 2.0, THESIS_RISE)
            (briggs < 0.99cap(BRIGGS_RISE) && thesis < 0.99cap(THESIS_RISE)) || continue
            checked += 1
            @test briggs ≈ two_thirds rtol = 6e-3
            @test thesis / two_thirds ≈ (6 / 1.6^3)^(1 / 3) rtol = 2e-3
        end
        @test checked ≥ 6          # the sweep must actually exercise the branch
        @test BRIGGS_RISE.combined_buoyancy == 2 * 0.6^2       # Briggs, β = 0.6
        @test 3 / 1.6^3 ≈ 0.7324 rtol = 1e-4                   # the same from the law
        @test THESIS_RISE.combined_buoyancy == 0.5
    end

    # The neutral momentum rise is 3 w₀D/u in Briggs (1969) Eq. 5.2, as
    # implemented by EPA ISC3 Eq. (1-16). The 2021 code had 1.5, unsourced.
    @testset "neutral momentum rise is Briggs" begin
        for w₀ in (5.0, 15.0), D in (1.0, 3.0), u in (2.0, 8.0)
            @test final_momentum_rise(0.0, w₀, D, u, -1e-6, BRIGGS_RISE) ≈ 3 * w₀ * D / u
            @test final_momentum_rise(0.0, w₀, D, u, -1e-6, THESIS_RISE) ≈ 1.5 * w₀ * D / u
        end
        @test BRIGGS_RISE.neutral_momentum == 3.0
    end

    # The stable final rise coefficient is 2.6 in Briggs and the Handbook on
    # Atmospheric Diffusion, and 2.4 in NRC XOQDOQ.
    @testset "stable rise coefficient, Briggs against XOQDOQ" begin
        @test BRIGGS_RISE.stable_final == 2.6
        @test XOQDOQ_RISE.stable_final == 2.4
        for F in (5.0, 200.0), u in (3.0, 8.0), S in (1e-4, 1e-3)
            b = final_buoyant_rise(F, u, S, BRIGGS_RISE)
            x = final_buoyant_rise(F, u, S, XOQDOQ_RISE)
            @test x ≤ b
        end
    end

    # Briggs (1973) open-country lateral dispersion, as tabulated in Hanna,
    # Briggs and Hosker, Handbook on Atmospheric Diffusion, DOE/TIC-11223
    # (1982), Table 4.5, attributing ATDL Contribution No. 79:
    # σ_y = a x (1 + 10⁻⁴x)^(−1/2) with a running
    # 0.22, 0.16, 0.11, 0.08, 0.06, 0.04 from class A to F. The table states
    # its own validity band as 10² < x < 10⁴ m.
    @testset "σ_y is Briggs 1973 open country" begin
        briggs_a = (0.22, 0.16, 0.11, 0.08, 0.06, 0.04)
        for (i, class) in enumerate(PASQUILL_CLASSES)
            @test lateral_coefficient(class) == briggs_a[i]
            for x in (10.0, 100.0, 1000.0, 10_000.0, 50_000.0)
                @test lateral_dispersion(x, class) ≈ briggs_a[i] * x * (1 + 1e-4 * x)^(-0.5) rtol =
                    1e-12
            end
        end
    end

    # Briggs distance to final rise: x_f = 14F^(5/8) below 55 m⁴/s³,
    # 34F^(2/5) above.
    @testset "distance to final rise is Briggs" begin
        for F in (1.0, 10.0, 54.9)
            @test buoyancy_transition_distance(F) ≈ 14 * F^(5 / 8) rtol = 1e-12
        end
        for F in (55.0, 200.0, 5000.0)
            @test buoyancy_transition_distance(F) ≈ 34 * F^(2 / 5) rtol = 1e-12
        end
    end

    # The neutral final buoyant rise is written here as
    # 1.6F^(1/3)(3.5x_f)^(2/3)/u. Briggs publishes it as 21.4F^(3/4)/u and
    # 38.7F^(3/5)/u. Substituting x_f shows these are one expression:
    # 1.6·49^(2/3) = 21.425 and 1.6·119^(2/3) = 38.71, which the literature
    # rounds. The agreement is therefore exact up to that rounding, and the
    # test asserts both the algebra and the numbers.
    @testset "neutral final rise is the published Briggs form" begin
        @test 1.6 * 49^(2 / 3) ≈ 21.425 rtol = 1e-4
        @test 1.6 * 119^(2 / 3) ≈ 38.71 rtol = 1e-4
        for F in (5.0, 20.0, 54.0), u in (3.0, 8.0)
            @test final_buoyant_rise(F, u, -1e-6) ≈ 21.4 * F^0.75 / u rtol = 2e-3
        end
        for F in (56.0, 200.0, 1000.0), u in (3.0, 8.0)
            @test final_buoyant_rise(F, u, -1e-6) ≈ 38.7 * F^0.6 / u rtol = 2e-3
        end
    end

    # Briggs stable final rise, 2.6[F/(u s)]^(1/3), and the two-thirds law
    # for the transitional rise.
    @testset "stable final rise and the two-thirds law" begin
        for F in (10.0, 100.0), u in (2.0, 6.0), S in (1e-4, 1e-3)
            stable = 2.6 * (F / (u * S))^(1 / 3)
            @test final_buoyant_rise(F, u, S) <= stable * (1 + 1e-12)
            # where the stable branch is the binding one, it is exact
            if stable < 1.6 * F^(1/3) * (3.5 * buoyancy_transition_distance(F))^(2/3) / u &&
               stable < 5.0 * F^(1/4) * S^(-3/8)
                @test final_buoyant_rise(F, u, S) ≈ stable rtol = 1e-12
            end
        end
        for F in (10.0, 100.0), u in (2.0, 6.0), x in (10.0, 100.0)
            two_thirds = 1.6 * F^(1 / 3) * x^(2 / 3) / u
            r = buoyant_rise(x, F, u, 1e-3)
            @test r ≈ min(two_thirds, final_buoyant_rise(F, u, 1e-3)) rtol = 1e-12
        end
    end

    # σ_z is *not* Briggs — it is the g(x)F(x) form of the normative, an
    # NRPB-R91-style scheme with an explicit roughness correction. It is not
    # asserted equal to Briggs; it is asserted to agree within the factor
    # that separates published σ schemes, with the expected sign: less
    # vertical spread than Briggs in unstable air, more in stable.
    @testset "σ_z brackets Briggs open country" begin
        briggs_z = Dict(
            PASQUILL_A => (x -> 0.20x),
            PASQUILL_B => (x -> 0.12x),
            PASQUILL_C => (x -> 0.08x * (1 + 0.0002x)^(-0.5)),
            PASQUILL_D => (x -> 0.06x * (1 + 0.0015x)^(-0.5)),
            PASQUILL_E => (x -> 0.03x * (1 + 0.0003x)^(-1)),
            PASQUILL_F => (x -> 0.016x * (1 + 0.0003x)^(-1)),
        )
        for class in PASQUILL_CLASSES, x in (100.0, 300.0, 1000.0, 3000.0, 10_000.0)
            ratio = vertical_dispersion(x, class, ROUGHNESS_PASTURE) / briggs_z[class](x)
            @test 0.4 < ratio < 1.6
        end
        # the sign of the difference, which is systematic
        @test vertical_dispersion(1000.0, PASQUILL_A, ROUGHNESS_PASTURE) <
              briggs_z[PASQUILL_A](1000.0)
        @test vertical_dispersion(1000.0, PASQUILL_F, ROUGHNESS_PASTURE) >
              briggs_z[PASQUILL_F](1000.0)
    end
end
