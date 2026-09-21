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
        indices = [sector_of(g, deg2rad(22.5 * (i - 1))) for i in 1:16]

        @test indices == collect(1:16)
        @test length(unique(indices)) == 16
        @test [sector_name(g, k) for k in indices] == collect(CARDINAL_16)

        for k in 1:16
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
        for k in 1:16
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
        for k in 1:16, offset in (-0.4, -0.1, 0.0, 0.1, 0.4)

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
        for j in 0:(m - 1)
            counts[sector_of(g, 2π * (j + 0.5) / m)] += 1
        end
        @test all(==(64), counts)
    end

    @testset "opposite" begin
        for n in (4, 8, 16, 36)
            g = SectorGrid(n)
            for k in 1:n
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
