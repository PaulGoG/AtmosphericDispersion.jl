using AtmosphericDispersion
using Test
using Aqua
using ExplicitImports
using JET

@testset "AtmosphericDispersion.jl" begin
    @testset "Static quality assurance" begin
        Aqua.test_all(AtmosphericDispersion)
        @test check_no_implicit_imports(AtmosphericDispersion) === nothing
        @test check_no_stale_explicit_imports(AtmosphericDispersion) === nothing
        @test check_all_qualified_accesses_via_owners(AtmosphericDispersion) === nothing
        JET.test_package(AtmosphericDispersion)
    end

    @testset "SectorGrid" begin
        @testset "construction" begin
            @test nsectors(SectorGrid(16)) == 16
            @test sector_width(SectorGrid(16)) ≈ deg2rad(22.5)
            @test sector_width(SectorGrid(8)) ≈ deg2rad(45.0)
            @test_throws ArgumentError SectorGrid(2)
            @test_throws ArgumentError SectorGrid(15)
        end

        @testset "bearing" begin
            @test bearing(0, 1) ≈ 0                # north
            @test bearing(1, 0) ≈ π / 2            # east
            @test bearing(0, -1) ≈ π               # south
            @test bearing(-1, 0) ≈ 3π / 2          # west
            @test bearing(1, 1) ≈ π / 4            # north-east
            @test bearing(0, 0) == 0
            @test 0 ≤ bearing(-1, -1) < 2π
        end

        # The 2021 implementation binned by truncating a ratio, which put every
        # cardinal direction on a sector edge and let one ulp of round-off decide
        # the bin: ESE and SE collapsed into one sector, SSW and SW into another,
        # and two sectors became unreachable. Each cardinal direction must land in
        # its own sector, at that sector's centre.
        @testset "the sixteen cardinal directions" begin
            g = SectorGrid(16)
            indices = [sector_of(g, deg2rad(22.5 * (i - 1))) for i = 1:16]

            @test indices == collect(1:16)
            @test length(unique(indices)) == 16
            @test [sector_name(g, k) for k in indices] == collect(CARDINAL_16)

            for k = 1:16
                @test sector_bearing(g, k) ≈ deg2rad(22.5 * (k - 1))
                @test sector_of(g, sector_bearing(g, k)) == k
            end
        end

        @testset "cardinal directions from Cartesian displacements" begin
            g = SectorGrid(16)
            # x east, y north; a receptor due east of the source is in sector 5.
            @test sector_of(g, 0.0, 1.0) == 1     # N
            @test sector_of(g, 1.0, 1.0) == 3     # NE
            @test sector_of(g, 1.0, 0.0) == 5     # E
            @test sector_of(g, 1.0, -1.0) == 7    # SE
            @test sector_of(g, 0.0, -1.0) == 9    # S
            @test sector_of(g, -1.0, -1.0) == 11  # SW
            @test sector_of(g, -1.0, 0.0) == 13   # W
            @test sector_of(g, -1.0, 1.0) == 15   # NW
            # Distance must not affect the sector.
            @test sector_of(g, 1e-9, 0.0) == sector_of(g, 1e9, 0.0)
        end

        @testset "sector boundaries" begin
            g = SectorGrid(16)
            half = sector_width(g) / 2
            for k = 1:16
                c = sector_bearing(g, k)
                lo, hi = sector_bounds(g, k)
                @test mod2pi(c - half) ≈ lo
                @test mod2pi(c + half) ≈ hi
                # Interior points either side of centre stay in the sector.
                @test sector_of(g, mod2pi(c - 0.99half)) == k
                @test sector_of(g, mod2pi(c + 0.99half)) == k
                # A bearing exactly on a boundary belongs to one of the two
                # sectors meeting there. The tie rule takes the higher, but a
                # boundary bearing is generally not representable, so round-off
                # may place it either side. Both answers are defensible.
                @test sector_of(g, mod2pi(c + half)) in (k, mod(k, 16) + 1)
            end
        end

        @testset "wrapping" begin
            g = SectorGrid(16)
            # Away from the boundaries the sector is invariant under whole turns.
            for k = 1:16, offset in (-0.4, -0.1, 0.0, 0.1, 0.4)
                β = sector_bearing(g, k) + offset * sector_width(g)
                @test sector_of(g, β) == k
                @test sector_of(g, β + 2π) == k
                @test sector_of(g, β - 2π) == k
                @test sector_of(g, β + 6π) == k
            end
            # Any bearing whatsoever lands in a valid sector.
            for β in range(-4π, 4π; length = 401)
                @test 1 ≤ sector_of(g, β) ≤ 16
            end
        end

        @testset "the sectors partition the circle" begin
            g = SectorGrid(16)
            # Sample midway between boundaries so no sample sits on one.
            counts = zeros(Int, 16)
            m = 16 * 64
            for j = 0:(m-1)
                counts[sector_of(g, 2π * (j + 0.5) / m)] += 1
            end
            @test all(==(64), counts)
        end

        @testset "opposite" begin
            for n in (4, 8, 16, 36)
                g = SectorGrid(n)
                for k = 1:n
                    o = opposite(g, k)
                    @test 1 ≤ o ≤ n
                    @test opposite(g, o) == k
                    @test o != k
                    # Opposite sectors are half a turn apart.
                    δ = abs(sector_bearing(g, k) - sector_bearing(g, o))
                    @test min(δ, 2π - δ) ≈ π
                end
            end
            g = SectorGrid(16)
            @test sector_name(g, opposite(g, 1)) == "S"
            @test sector_name(g, opposite(g, 5)) == "W"
            @test_throws BoundsError opposite(g, 0)
            @test_throws BoundsError opposite(g, 17)
        end

        @testset "naming" begin
            g = SectorGrid(16)
            @test sector_name(g, 1) == "N"
            @test sector_name(g, 5) == "E"
            @test sector_name(g, 9) == "S"
            @test sector_name(g, 13) == "W"
            @test sector_name(SectorGrid(8), 3) == "90.0°"
        end
    end

    @testset "PasquillClass" begin
        @test length(PASQUILL_CLASSES) == 6
        @test PASQUILL_CLASSES[1] === PASQUILL_A
        @test PASQUILL_CLASSES[6] === PASQUILL_F
        @test pasquill('D') === PASQUILL_D
        @test pasquill('d') === PASQUILL_D
        @test pasquill("F") === PASQUILL_F
        @test pasquill(PASQUILL_B) === PASQUILL_B
        @test_throws ArgumentError pasquill('G')
        @test_throws ArgumentError pasquill("DE")
        for (i, c) in enumerate(PASQUILL_CLASSES)
            @test classindex(c) == i
            @test pasquill(letter(c)) === c
        end
    end

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
            @test_throws ArgumentError WindRose(
                g,
                f,
                BlowingFrom();
                stability = fill(1 / 3, 6),
            )
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
end
