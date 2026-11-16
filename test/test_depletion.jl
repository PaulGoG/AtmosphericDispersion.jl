@testset "Depletion" begin
    site = REFERENCE_SITE
    inert = INERT

    @testset "nuclides" begin
        @test TRITIATED_WATER.decay_constant == TRITIUM_DECAY_CONSTANT
        @test TRITIUM_GAS.decay_constant == TRITIUM_DECAY_CONSTANT
        # HTO deposits an order of magnitude more readily than HT: exactly
        # a factor of ten in the low velocity, sixteen in the high.
        @test TRITIATED_WATER.deposition_velocity.low ≈
              10 * TRITIUM_GAS.deposition_velocity.low
        @test TRITIATED_WATER.deposition_velocity.high ≈
              16 * TRITIUM_GAS.deposition_velocity.high
        # Tritium's half-life is about 12.3 years.
        @test log(2) / TRITIUM_DECAY_CONSTANT / (365.25 * 86400) ≈ 12.3 atol = 0.05
        @test_throws ArgumentError DepositionVelocity(0.5, 0.1)
        @test_throws ArgumentError DepositionVelocity(-0.1, 0.1)
        @test_throws ArgumentError Nuclide(;
            name = "x",
            decay_constant = -1.0,
            deposition_velocity = DepositionVelocity(0.0, 0.0),
        )
        @test occursin("HTO", sprint(show, TRITIATED_WATER))
    end

    @testset "decay" begin
        u = transport_wind_speed(site, PASQUILL_D)
        @test decay_factor(0.0, u, TRITIATED_WATER) == 1
        @test decay_factor(1000.0, u, TRITIATED_WATER) ==
              exp(-TRITIUM_DECAY_CONSTANT * 1000 / u)
        @test decay_factor(1000.0, u, inert) == 1
        # Tritium is far too long-lived for transit decay to matter.
        @test decay_factor(1e5, u, TRITIATED_WATER) > 0.9999
        @test issorted(
            [decay_factor(x, u, TRITIATED_WATER) for x in (0.0, 1e4, 1e6, 1e8)];
            rev = true,
        )
        @test_throws DomainError decay_factor(-1.0, u, TRITIATED_WATER)
        @test_throws DomainError decay_factor(100.0, 0.0, TRITIATED_WATER)
    end

    @testset "dry depletion" begin
        x_floor = AtmosphericDispersion.DEPLETION_INTEGRAL_FLOOR
        for class in PASQUILL_CLASSES
            @test depletion_integral(x_floor, site, class) == 0
            @test depletion_integral(0.5, site, class) == 0
            # The integral accumulates with distance.
            I = [depletion_integral(x, site, class) for x in (10.0, 1e3, 1e4, 1e5)]
            @test issorted(I)
            @test all(≥(0), I)
            # A surviving fraction, so in (0, 1] and falling with distance.
            f = [dry_depletion_factor(x, site, class, TRITIATED_WATER)
                 for
                 x in (10.0, 1e3, 1e4, 1e5)]
            @test all(v -> 0 < v ≤ 1, f)
            @test issorted(f; rev = true)
            # No deposition velocity, no depletion.
            @test dry_depletion_factor(1e4, site, class, inert) == 1
            # HT deposits less than HTO, so more of it survives.
            @test dry_depletion_factor(1e4, site, class, TRITIUM_GAS) >
                  dry_depletion_factor(1e4, site, class, TRITIATED_WATER)
        end
    end

    @testset "WashoutEvent" begin
        event = WashoutEvent(; duration = 3600)
        @test event.precipitation === PRECIPITATION_RAIN
        @test event.rate == first(PRECIPITATION_RATES)
        @test event.duration == 3600.0
        @test event.model === WASHOUT_NORMATIVE
        # The washout table is defined only at the tabulated intensities.
        @test_throws ArgumentError WashoutEvent(; duration = 3600.0, rate = 2.0)
        @test_throws ArgumentError WashoutEvent(; duration = 3600.0, rate = 0.0)
        @test_throws ArgumentError WashoutEvent(; duration = -1.0)
        @test_throws ArgumentError WashoutEvent(;
            duration = -1.0,
            precipitation = PRECIPITATION_RAIN,
            rate = 1.0,
        )
        # The event looks its own coefficients up.
        for p in PRECIPITATION_TYPES,
            r in PRECIPITATION_RATES,
            m in (WASHOUT_NORMATIVE, WASHOUT_HTO)
            e = WashoutEvent(; duration = 600.0, precipitation = p, rate = r, model = m)
            for species in (WASHOUT_TRITIUM_IODINE, WASHOUT_OTHER_NUCLIDES)
                @test washout_coefficients(e, species) ==
                      washout_coefficients(p, r, m, species)
            end
            @test washout_coefficients(e) == washout_coefficients(p, r, m)
        end
    end

    @testset "wet depletion" begin
        wet(t, p, r) = wet_depletion_factor(WashoutEvent(; duration = t, precipitation = p, rate = r))
        @test wet(0.0, PRECIPITATION_RAIN, 1.0) == 1
        for p in PRECIPITATION_TYPES, r in PRECIPITATION_RATES

            f = [wet(t, p, r) for t in (0.0, 600.0, 3600.0, 86400.0)]
            @test all(v -> 0 < v ≤ 1, f)
            @test issorted(f; rev = true)
        end
        # Heavier rain washes out faster; snow barely at all.
        @test wet(3600.0, PRECIPITATION_RAIN, 5.0) < wet(3600.0, PRECIPITATION_RAIN, 0.5)
        @test wet(3600.0, PRECIPITATION_SNOW, 5.0) > wet(3600.0, PRECIPITATION_RAIN, 0.5)
    end

    # The composition is multiplicative, not additive. With nothing removed
    # every factor is one and so is the product. The 2021 code added the wet
    # and dry factors and returned two in this limit.
    @testset "the no-depletion limit is one" begin
        no_rain = WashoutEvent(; duration = 0.0)
        for class in PASQUILL_CLASSES, x in (10.0, 1e3, 1e5)

            @test depletion_factor(x, site, class, inert) == 1
            # Rain that lasts no time removes nothing either.
            @test depletion_factor(x, site, class, inert; washout = no_rain) == 1
        end
    end

    @testset "washout acts for the travel time, on the share of the release in rain" begin
        u, x = 6.0, 6000.0
        Λ = 1e-5        # NSR-23 Table 7, tritium and iodine, rain at 1 mm/h, low
        whole = WashoutEvent(; duration = 3600.0, rate = 1.0)
        quarter = WashoutEvent(; duration = 900.0, release_duration = 3600.0, rate = 1.0)
        dry = WashoutEvent(; duration = 0.0, rate = 1.0)
        @test precipitation_fraction(whole) == 1
        @test precipitation_fraction(quarter) == 0.25
        @test precipitation_fraction(dry) == 0
        @test wet_depletion_factor(x, u, whole) ≈ exp(-Λ * x / u)
        @test wet_depletion_factor(x, u, dry) == 1
        # The quarter-hour of rain cannot scavenge for the 1000 s of travel:
        # what was released into it is under rain for 900 s at most.
        @test wet_depletion_factor(x, u, quarter) ≈ 0.75 + 0.25 * exp(-Λ * 900)
        @test wet_depletion_factor(4800.0, u, quarter) ≈ 0.75 + 0.25 * exp(-Λ * 800)
        # Continuous at x = u t_w, and constant beyond.
        @test wet_depletion_factor(5400.0, u, quarter) ≈
              wet_depletion_factor(5400.0 + 1e-6, u, quarter)
        @test wet_depletion_factor(5e4, u, quarter) ==
              wet_depletion_factor(6000.0, u, quarter)
        # Far out, precipitation through the whole release gives the 2021 form.
        @test wet_depletion_factor(1e7, u, whole) ≈ wet_depletion_factor(whole)
        # At the source nothing has been scavenged yet, however long it rains.
        @test wet_depletion_factor(0.0, u, whole) == 1
        # Rain longer than the release covers all of it.
        @test precipitation_fraction(
            WashoutEvent(; duration = 7200.0, release_duration = 3600.0),) == 1
        @test_throws ArgumentError WashoutEvent(; duration = 10.0, release_duration = 0.0)
        @test_throws ArgumentError WashoutEvent(; duration = 10.0, release_duration = -1.0)
    end

    @testset "composition" begin
        for class in PASQUILL_CLASSES
            x, t = 1e4, 3600.0
            rain = WashoutEvent(; duration = t, precipitation = PRECIPITATION_RAIN, rate = 1.0)
            u = transport_wind_speed(site, class)
            # Rain through the whole release: the material is scavenged for
            # its travel time x/u with the low coefficient of NSR-23 Table 7,
            # 10⁻⁵ s⁻¹ at 1 mm/h.
            expected = decay_factor(x, u, TRITIATED_WATER) *
                       dry_depletion_factor(x, site, class, TRITIATED_WATER) *
                       exp(-1e-5 * x / u)
            @test depletion_factor(x, site, class, TRITIATED_WATER; washout = rain) ≈
                  expected
            @test 0 < expected ≤ 1
            # Rain can only remove more.
            @test depletion_factor(
                x,
                site,
                class,
                TRITIATED_WATER;
                washout = WashoutEvent(; duration = t),
            ) < depletion_factor(x, site, class, TRITIATED_WATER)
        end
    end

    @testset "depletion inside the long-term class sum" begin
        g = SectorGrid(16)
        rose = WindRose(
            g,
            fill(1 / 16, 16),
            BlowingToward();
            stability = normalised_stability(),
        )
        bare = dilution_long_term(0.0, -5000.0, site, rose)
        # An inert species leaves the field untouched: the depletion factor
        # is exactly one, not six as the unweighted sum of the 2021 code
        # would have given.
        @test dilution_long_term(0.0, -5000.0, site, rose; nuclide = inert) == bare
        depleted = dilution_long_term(0.0, -5000.0, site, rose; nuclide = TRITIATED_WATER)
        @test 0 < depleted < bare
        @test depleted / bare > 0.9        # tritium depletes slowly
        # Washout removes more still.
        @test dilution_long_term(
            0.0,
            -5000.0,
            site,
            rose;
            nuclide = TRITIATED_WATER,
            washout = WashoutEvent(; duration = 3600.0),
        ) < depleted
    end
end
