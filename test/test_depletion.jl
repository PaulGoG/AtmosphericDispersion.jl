@testset "Depletion" begin
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
    site = Site(; source = stack, atmosphere = air)
    inert = Nuclide(;
        name = "inert",
        decay_constant = 0.0,
        deposition_velocity = DepositionVelocity(0.0, 0.0),
    )

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
        for class in PASQUILL_CLASSES
            @test depletion_integral(DEPLETION_INTEGRAL_FLOOR, site, class) == 0
            @test depletion_integral(0.5, site, class) == 0
            # The integral accumulates with distance.
            I = [depletion_integral(x, site, class) for x in (10.0, 1e3, 1e4, 1e5)]
            @test issorted(I)
            @test all(≥(0), I)
            # A surviving fraction, so in (0, 1] and falling with distance.
            f = [
                dry_depletion_factor(x, site, class, TRITIATED_WATER) for
                x in (10.0, 1e3, 1e4, 1e5)
            ]
            @test all(v -> 0 < v ≤ 1, f)
            @test issorted(f; rev = true)
            # No deposition velocity, no depletion.
            @test dry_depletion_factor(1e4, site, class, inert) == 1
            # HT deposits less than HTO, so more of it survives.
            @test dry_depletion_factor(1e4, site, class, TRITIUM_GAS) >
                  dry_depletion_factor(1e4, site, class, TRITIATED_WATER)
        end
    end

    @testset "wet depletion" begin
        @test wet_depletion_factor(0.0, PRECIPITATION_RAIN, 1.0) == 1
        for p in PRECIPITATION_TYPES, r in PRECIPITATION_RATES
            f = [wet_depletion_factor(t, p, r) for t in (0.0, 600.0, 3600.0, 86400.0)]
            @test all(v -> 0 < v ≤ 1, f)
            @test issorted(f; rev = true)
        end
        # Heavier rain washes out faster; snow barely at all.
        @test wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 5.0) <
              wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 0.5)
        @test wet_depletion_factor(3600.0, PRECIPITATION_SNOW, 5.0) >
              wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 0.5)
        @test_throws DomainError wet_depletion_factor(-1.0, PRECIPITATION_RAIN, 1.0)
    end

    # The composition is multiplicative, not additive. With nothing removed
    # at all every factor is one and so is the product. The 2021 code added
    # the wet and dry factors and returned two in this limit.
    @testset "the no-depletion limit is one" begin
        for class in PASQUILL_CLASSES, x in (10.0, 1e3, 1e5)
            @test depletion_factor(x, site, class, inert) == 1
            @test depletion_factor(x, site, class, inert; washout_duration = 0.0) == 1
        end
    end

    @testset "composition" begin
        for class in PASQUILL_CLASSES
            x, t = 1e4, 3600.0
            u = transport_wind_speed(site, class)
            expected =
                decay_factor(x, u, TRITIATED_WATER) *
                dry_depletion_factor(x, site, class, TRITIATED_WATER) *
                wet_depletion_factor(t, PRECIPITATION_RAIN, 1.0)
            @test depletion_factor(
                x,
                site,
                class,
                TRITIATED_WATER;
                washout_duration = t,
                precipitation = PRECIPITATION_RAIN,
                rate = 1.0,
            ) ≈ expected
            @test 0 < expected ≤ 1
            # Rain can only remove more.
            @test depletion_factor(x, site, class, TRITIATED_WATER; washout_duration = t) <
                  depletion_factor(x, site, class, TRITIATED_WATER)
        end
    end

    @testset "depletion inside the long-term class sum" begin
        g = SectorGrid(16)
        stab = [0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158]
        stab ./= sum(stab)
        rose = WindRose(g, fill(1 / 16, 16), BlowingToward(); stability = stab)
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
            washout_duration = 3600.0,
        ) < depleted
    end
end
