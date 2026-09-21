@testset "WindRose" begin
    g = SectorGrid(16)
    uniform_stability = fill(1 / 6, 6)

    @testset "the convention is a half-turn rotation" begin
        f = zeros(16)
        f[1] = 1.0   # all the wind is in sector 1 (north)
        from = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
        toward = WindRose(g, f, BlowingToward(); stability = uniform_stability)

        # Wind from the north transports towards the south, sector 9.
        @test frequency_toward(from, 9) == 1.0
        @test frequency_toward(from, 1) == 0.0
        @test frequency_from(from, 1) == 1.0

        # Declared the other way, the same numbers mean the opposite field.
        @test frequency_toward(toward, 1) == 1.0
        @test frequency_toward(toward, 9) == 0.0
        @test frequency_from(toward, 9) == 1.0

        # The two readings of one table are exactly opposite everywhere.
        for k = 1:16
            @test frequency_toward(from, k) == frequency_toward(toward, opposite(g, k))
            @test frequency_from(from, k) == frequency_toward(from, opposite(g, k))
        end
    end

    @testset "a symmetric rose is invariant under the convention" begin
        f = fill(1 / 16, 16)
        a = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
        b = WindRose(g, f, BlowingToward(); stability = uniform_stability)
        @test all(k -> frequency_toward(a, k) == frequency_toward(b, k), 1:16)
    end

    @testset "stability fractions" begin
        f = fill(1 / 16, 16)
        s = [0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158]  # the 2021 input table
        rose = WindRose(g, f, BlowingFrom(); stability = s)
        for k = 1:16
            @test sum(stability_fraction(rose, k, c) for c in PASQUILL_CLASSES) ≈ 1 atol =
                1e-3
            @test stability_fraction(rose, k, PASQUILL_D) == 0.488
        end

        # A per-sector matrix must follow the sector rotation of its rose.
        m = repeat(reshape(s, 1, 6), 16, 1)
        m[1, :] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        rose_m = WindRose(g, f, BlowingToward(); stability = m)
        @test stability_fraction(rose_m, 1, PASQUILL_D) == 1.0
        @test stability_fraction(rose_m, 2, PASQUILL_D) == 0.488
    end

    @testset "validation" begin
        f = fill(1 / 16, 16)
        @test_throws DimensionMismatch WindRose(
            g,
            fill(1 / 8, 8),
            BlowingFrom();
            stability = uniform_stability,
        )
        @test_throws DimensionMismatch WindRose(
            g,
            f,
            BlowingFrom();
            stability = fill(1 / 5, 5),
        )
        @test_throws DimensionMismatch WindRose(
            g,
            f,
            BlowingFrom();
            stability = fill(1 / 6, 8, 6),
        )
        # Raw counts rather than a normalised distribution.
        @test_throws ArgumentError WindRose(
            g,
            fill(10.0, 16),
            BlowingFrom();
            stability = uniform_stability,
        )
        bad = copy(f)
        bad[1] = -bad[2]
        bad[2] *= 2
        @test_throws ArgumentError WindRose(
            g,
            bad,
            BlowingFrom();
            stability = uniform_stability,
        )
        nan = copy(f)
        nan[1] = NaN
        @test_throws ArgumentError WindRose(
            g,
            nan,
            BlowingFrom();
            stability = uniform_stability,
        )
        @test_throws ArgumentError WindRose(g, f, BlowingFrom(); stability = fill(1 / 3, 6))
    end

    @testset "accessors" begin
        f = fill(1 / 16, 16)
        rose = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
        @test grid(rose) === g
        @test nsectors(rose) == 16
        @test_throws BoundsError frequency_toward(rose, 0)
        @test_throws BoundsError stability_fraction(rose, 17, PASQUILL_A)
        @test occursin("WindRose", sprint(show, MIME"text/plain"(), rose))
    end
end
