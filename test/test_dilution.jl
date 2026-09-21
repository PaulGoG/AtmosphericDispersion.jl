@testset "Dilution" begin
    site = REFERENCE_SITE

    @testset "plume frame" begin
        # Wind from the north transports to the south.
        x, y = plume_frame(0.0, -1000.0, 0.0)
        @test x ≈ 1000 && abs(y) < 1e-9
        @test plume_frame(0.0, 1000.0, 0.0)[1] ≈ -1000        # upwind
        # Wind from the west transports to the east.
        x, y = plume_frame(1000.0, 0.0, 3π / 2)
        @test x ≈ 1000 && abs(y) < 1e-9
        # The frame is orthonormal: it preserves distance.
        for β in range(0, 2π; length = 17), (e, n) in ((300.0, 400.0), (-120.0, 50.0))

            x, y = plume_frame(e, n, β)
            @test hypot(x, y) ≈ hypot(e, n)
        end
    end

    @testset "instantaneous field" begin
        # Nothing upwind.
        @test dilution_instantaneous(0.0, 1000.0, 0.0, site, PASQUILL_D, 0.0) == 0
        # Peaked on the axis and symmetric across it.
        onaxis = dilution_instantaneous(0.0, -1000.0, 0.0, site, PASQUILL_D, 0.0)
        @test onaxis > 0
        for off in (50.0, 200.0)
            left = dilution_instantaneous(off, -1000.0, 0.0, site, PASQUILL_D, 0.0)
            right = dilution_instantaneous(-off, -1000.0, 0.0, site, PASQUILL_D, 0.0)
            @test left ≈ right
            @test left < onaxis
        end
        # Falls off downwind, far from the source.
        far = [dilution_instantaneous(0.0, -r, 0.0, site, PASQUILL_D, 0.0)
               for
               r in (2_000.0, 5_000.0, 20_000.0)]
        @test issorted(far; rev = true)
        # Equivariant under a common rotation of receptor and wind.
        for β in range(0, 2π; length = 13)
            rotated = dilution_instantaneous(
                1000sin(β + π),
                1000cos(β + π),
                0.0,
                site,
                PASQUILL_D,
                β,
            )
            @test rotated ≈ onaxis
        end
        @test_throws DomainError dilution_instantaneous(
            0.0,
            -1000.0,
            -1.0,
            site,
            PASQUILL_D,
            0.0,
        )
    end

    @testset "extended field" begin
        # Uniform across the sector, zero outside it.
        onaxis = dilution_extended(0.0, -1000.0, site, PASQUILL_D, 0.0)
        @test onaxis > 0
        inside = dilution_extended(90.0, -1000.0, site, PASQUILL_D, 0.0)
        @test inside ≈ onaxis        # crosswind-uniform within the sector
        @test dilution_extended(1000.0, 0.0, site, PASQUILL_D, 0.0) == 0   # 90° off
        @test dilution_extended(0.0, 1000.0, site, PASQUILL_D, 0.0) == 0   # upwind
        # A wider sector spreads the same material further, so dilutes more.
        @test dilution_extended(0.0, -1000.0, site, PASQUILL_D, 0.0, SectorGrid(8)) < onaxis
    end

    @testset "long-term field" begin
        g = SectorGrid(16)
        stab = normalised_stability()

        # A uniform rose gives a rotationally symmetric field.
        uniform = WindRose(g, fill(1 / 16, 16), BlowingToward(); stability = stab)
        values = [dilution_long_term(
                      1000sin(sector_bearing(g, k)),
                      1000cos(sector_bearing(g, k)),
                      site,
                      uniform,
                  ) for k in 1:16]
        @test all(v -> v ≈ first(values), values)
        @test first(values) > 0

        # The two conventions give fields that are exact mirror images,
        # sector by sector. This is the defect, isolated.
        f = [
            0.07849,
            0.07849,
            0.07221,
            0.07378,
            0.07849,
            0.05965,
            0.05651,
            0.0471,
            0.05181,
            0.05024,
            0.04867,
            0.04553,
            0.05495,
            0.06436,
            0.06279,
            0.07692,
        ]
        f ./= sum(f)
        from = WindRose(g, f, BlowingFrom(); stability = stab)
        toward = WindRose(g, f, BlowingToward(); stability = stab)
        χ(rose, k) = dilution_long_term(
            1000sin(sector_bearing(g, k)),
            1000cos(sector_bearing(g, k)),
            site,
            rose,
        )
        for k in 1:16
            @test χ(from, k) ≈ χ(toward, opposite(g, k))
        end
        # And they disagree substantially where the rose is asymmetric.
        ratios = [χ(toward, k) / χ(from, k) for k in 1:16]
        @test maximum(ratios) > 1.4
        @test minimum(ratios) < 0.7
        # The most exposed sector is exactly opposite between the two.
        @test argmax([χ(toward, k) for k in 1:16]) ==
              opposite(g, argmax([χ(from, k) for k in 1:16]))

        @test dilution_long_term(0.0, 0.0, site, toward) == 0
        # A sector the wind never blows towards receives nothing.
        single = zeros(16)
        single[1] = 1.0
        only_north = WindRose(g, single, BlowingToward(); stability = stab)
        @test χ(only_north, 1) > 0
        @test χ(only_north, 9) == 0
    end

    @testset "long term is the rose-weighted extended regime" begin
        # With every hour spent in one class, the long-term factor on a
        # sector axis is the extended factor times the frequency of
        # transport towards that sector — identically, whatever the vertical
        # profile. It fails if the two regimes treat the mixing lid
        # differently, which they once did: the lid reached the extended
        # regime and the depletion integral but not the long-term sum.
        g = SectorGrid(16)
        F_k = fill(1 / 16, 16)
        lid_off = UNBOUNDED_SITE
        for (i, class) in enumerate(PASQUILL_CLASSES)
            stab = zeros(6)
            stab[i] = 1.0
            rose = WindRose(g, F_k, BlowingFrom(); stability = stab)
            for s in (site, lid_off), r in (2e3, 1e4, 2e4, 5e4), k in (1, 6, 12)
                β = sector_bearing(g, k)
                east, north = r * sin(β), r * cos(β)
                from = sector_bearing(g, opposite(g, k))
                @test dilution_long_term(east, north, s, rose) ≈
                      dilution_extended(east, north, s, class, from, g) / 16 rtol = 1e-12
            end
        end
        # The lid is felt where it should be: class A, whose σ_z reaches its
        # 1300 m depth within 20 km, and not class D within the same range.
        only_A = WindRose(g, F_k, BlowingFrom(); stability = [1.0, 0, 0, 0, 0, 0])
        only_D = WindRose(g, F_k, BlowingFrom(); stability = [0, 0, 0, 1.0, 0, 0])
        ratio(rose, r) = dilution_long_term(0.0, -r, site, rose) /
                         dilution_long_term(0.0, -r, lid_off, rose)
        @test ratio(only_A, 2e4) ≈ 1.44 atol = 0.01
        @test ratio(only_A, 5e4) ≈ 2.23 atol = 0.01
        @test ratio(only_D, 2e4) ≈ 1 atol = 1e-3
    end

    # Which row of the washout table applies is a property of the species, so
    # rain with nothing named for it to act on is refused rather than ignored.
    @testset "washout needs a nuclide" begin
        rain = WashoutEvent(; duration = 3600.0)
        rose = WindRose(
            SectorGrid(16),
            fill(1 / 16, 16),
            BlowingToward();
            stability = normalised_stability(),
        )
        @test_throws ArgumentError dilution_instantaneous(
            0.0,
            -1000.0,
            0.0,
            site,
            PASQUILL_D,
            0.0;
            washout = rain,
        )
        @test_throws ArgumentError dilution_extended(
            0.0,
            -1000.0,
            site,
            PASQUILL_D,
            0.0;
            washout = rain,
        )
        @test_throws ArgumentError dilution_long_term(
            0.0,
            -1000.0,
            site,
            rose;
            washout = rain,
        )
    end
end
