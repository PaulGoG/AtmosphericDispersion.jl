@testset "Deposition and resuspension" begin
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
    Q = 1e12

    @testset "dry deposition" begin
        @test dry_deposition(0.0, TRITIATED_WATER) == 0
        @test dry_deposition(1000.0, TRITIATED_WATER) ==
              TRITIATED_WATER.deposition_velocity.high * 1000
        # Deposition takes the high velocity, depletion the low one.
        @test dry_deposition(1000.0, TRITIATED_WATER) >
              TRITIATED_WATER.deposition_velocity.low * 1000
        # HT deposits far less than HTO.
        @test dry_deposition(1000.0, TRITIUM_GAS) < dry_deposition(1000.0, TRITIATED_WATER)
        @test_throws DomainError dry_deposition(-1.0, TRITIATED_WATER)
    end

    @testset "wet deposition" begin
        ω = wet_deposition(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
        )
        @test ω > 0
        # Nothing upwind, nothing for a null release.
        @test wet_deposition(
            0.0,
            1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
        ) == 0
        @test wet_deposition(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = 0.0,
        ) == 0
        # Linear in released activity.
        @test wet_deposition(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = 2Q,
            rate = 1.0,
        ) ≈ 2ω
        # Peaked on the axis, symmetric across it.
        off = wet_deposition(
            300.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
        )
        @test off < ω
        @test off ≈ wet_deposition(
            -300.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
        )
        # Heavier rain deposits more; snow far less.
        @test wet_deposition(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
            rate = 5.0,
        ) > ω
        @test wet_deposition(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
            precipitation = PRECIPITATION_SNOW,
        ) < ω / 100

        # The column normalisation is √(2π) Σ_y u, which is what integrating
        # the Gaussian plume over all z leaves. The 2021 code wrote √2 π,
        # larger by exactly √π, understating deposition by that factor.
        Λ = washout_coefficients(PRECIPITATION_RAIN, 1.0).high
        u = transport_wind_speed(site, PASQUILL_D)
        x, _ = plume_frame(0.0, -1000.0, 0.0)
        Σy = corrected_lateral_dispersion(x, site, PASQUILL_D)
        expected = Λ * Q * decay_factor(x, u, TRITIATED_WATER) / (sqrt(2π) * Σy * u)
        @test ω ≈ expected
        @test sqrt(2) * π / sqrt(2π) ≈ sqrt(π)
        @test ω ≈
              sqrt(π) * Λ * Q * decay_factor(x, u, TRITIATED_WATER) / (sqrt(2) * π * Σy * u)
    end

    @testset "sector-averaged wet deposition" begin
        ω = wet_deposition_sector(
            1000.0,
            site,
            PASQUILL_D,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
        )
        @test ω > 0
        @test wet_deposition_sector(0.0, site, PASQUILL_D, TRITIATED_WATER; activity = Q) ==
              0
        # Spread over an arc that grows with distance, so falls as 1/r.
        @test wet_deposition_sector(
            2000.0,
            site,
            PASQUILL_D,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
        ) ≈ ω / 2 rtol = 1e-6
        # A wider sector spreads the same material further.
        @test wet_deposition_sector(
            1000.0,
            site,
            PASQUILL_D,
            TRITIATED_WATER;
            activity = Q,
            rate = 1.0,
            sectors = SectorGrid(8),
        ) ≈ ω / 2 rtol = 1e-6
    end

    @testset "resuspension" begin
        c = RESUSPENSION_COEFFICIENTS
        @test resuspension_factor(0.0) ≈ c.A + c.B
        # Falls monotonically, and stays positive.
        ks = [resuspension_factor(t) for t in (0.0, 1.0, 30.0, 365.0, 3650.0)]
        @test issorted(ks; rev = true)
        @test all(>(0), ks)
        # About a factor 38 over the first year, four orders over ten.
        @test resuspension_factor(0.0) / resuspension_factor(365.0) ≈ 38 atol = 1
        @test resuspension_factor(0.0) / resuspension_factor(3650.0) > 1e4
        # The slow term is all that survives in the long run.
        @test resuspension_factor(1e5) ≈ c.B * exp(-c.λ₂ * 1e5)
        @test_throws DomainError resuspension_factor(-1.0)

        @test resuspended_concentration(0.0, 10.0) == 0
        @test resuspended_concentration(1000.0, 0.0) ≈ 1000 * (c.A + c.B)
        @test_throws DomainError resuspended_concentration(-1.0, 10.0)
    end

    @testset "depletion reaches the short-range regimes too" begin
        inert = Nuclide(;
            name = "inert",
            decay_constant = 0.0,
            deposition_velocity = DepositionVelocity(0.0, 0.0),
        )
        bare = dilution_instantaneous(0.0, -5000.0, 0.0, site, PASQUILL_D, 0.0)
        @test dilution_instantaneous(
            0.0,
            -5000.0,
            0.0,
            site,
            PASQUILL_D,
            0.0;
            nuclide = inert,
        ) == bare
        @test 0 <
              dilution_instantaneous(
                  0.0,
                  -5000.0,
                  0.0,
                  site,
                  PASQUILL_D,
                  0.0;
                  nuclide = TRITIATED_WATER,
              ) <
              bare
        bare_e = dilution_extended(0.0, -5000.0, site, PASQUILL_D, 0.0)
        @test dilution_extended(0.0, -5000.0, site, PASQUILL_D, 0.0; nuclide = inert) ==
              bare_e
        @test 0 <
              dilution_extended(
                  0.0,
                  -5000.0,
                  site,
                  PASQUILL_D,
                  0.0;
                  nuclide = TRITIATED_WATER,
              ) <
              bare_e
    end

    @testset "the washout model and species reach every kernel" begin
        # Reference values are read off NSR-23 Table 7 at 1 mm/h and Ogram
        # Eq. (38), not off the package: rain (low, high) is (1e-5, 2e-4) for
        # tritium and iodine and (2e-5, 3e-4) for every other nuclide; the
        # tritium snow column is (1e-7, 4e-7).
        other = Nuclide(;
            name = "other",
            decay_constant = TRITIUM_DECAY_CONSTANT,
            deposition_velocity = TRITIATED_WATER.deposition_velocity,
            washout_species = WASHOUT_OTHER_NUCLIDES,
        )
        @test TRITIATED_WATER.washout_species === WASHOUT_TRITIUM_IODINE
        T = 3600.0
        wash = (; washout_duration = T, precipitation = PRECIPITATION_RAIN, rate = 1.0)

        # Depletion takes the low coefficient, so the species moves every
        # dilution regime by exp(−ΔΛ_low T) and nothing else.
        expected = exp(-(2e-5 - 1e-5) * T)
        g = SectorGrid(16)
        # One class, so the depletion ratio is not blurred by the class sum.
        rose = WindRose(
            g,
            fill(1 / 16, 16),
            BlowingFrom();
            stability = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        )
        inst(n) = dilution_instantaneous(
            0.0,
            -5000.0,
            0.0,
            site,
            PASQUILL_D,
            0.0;
            nuclide = n,
            wash...,
        )
        ext(n) =
            dilution_extended(0.0, -5000.0, site, PASQUILL_D, 0.0; nuclide = n, wash...)
        lt(n) = dilution_long_term(0.0, -5000.0, site, rose; nuclide = n, wash...)
        for χ in (inst, ext, lt)
            @test χ(other) / χ(TRITIATED_WATER) ≈ expected rtol = 1e-12
        end

        # Deposition takes the high coefficient as well: 3e-4 against 2e-4.
        ω(n; kw...) =
            wet_deposition(0.0, -5000.0, site, PASQUILL_D, 0.0, n; activity = Q, kw...)
        ωs(n; kw...) =
            wet_deposition_sector(5000.0, site, PASQUILL_D, n; activity = Q, kw...)
        for f in (ω, ωs)
            @test f(other; wash...) / f(TRITIATED_WATER; wash...) ≈ 1.5 * expected rtol =
                1e-12
        end

        # The model changes snow only, and only the high coefficient: Ogram
        # at 1 mm/h is 1.2e-4 + 3.0e-4 against the tabulated 4e-7.
        snow = (; precipitation = PRECIPITATION_SNOW, rate = 1.0)
        for f in (ω, ωs)
            @test f(TRITIATED_WATER; snow..., washout_model = WASHOUT_HTO) /
                  f(TRITIATED_WATER; snow...) ≈ 4.2e-4 / 4e-7 rtol = 1e-12
            @test f(TRITIATED_WATER; wash..., washout_model = WASHOUT_HTO) ==
                  f(TRITIATED_WATER; wash...)
        end
        # The low coefficient is shared, so depletion does not see the model.
        @test inst(TRITIATED_WATER) == dilution_instantaneous(
            0.0,
            -5000.0,
            0.0,
            site,
            PASQUILL_D,
            0.0;
            nuclide = TRITIATED_WATER,
            wash...,
            washout_model = WASHOUT_HTO,
        )
    end
end
