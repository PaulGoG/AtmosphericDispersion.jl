@testset "Buildings" begin
    @testset "an empty envelope is the identity" begin
        e = BuildingEnvelope()
        @test equivalent_height(e) == 0
        @test equivalent_area(e) == 0
        for σ in (1.0, 50.0, 500.0), H in (0.0, 10.0, 100.0)

            @test wake_broadened(σ, H, e) == σ
        end
    end

    @testset "only buildings close enough count" begin
        near = Building(; east = 30.0, north = 0.0, height = 20.0, frontal_area = 400.0)
        far = Building(; east = 500.0, north = 0.0, height = 20.0, frontal_area = 400.0)
        @test equivalent_height(BuildingEnvelope([near])) == 20.0
        # 500 m away but only 20 m tall: outside 3 x its own height.
        @test equivalent_height(BuildingEnvelope([far])) == 0
        @test equivalent_height(BuildingEnvelope([near, far])) == 20.0
    end

    @testset "the tallest building in the vicinity governs" begin
        a = Building(; east = 10.0, north = 0.0, height = 30.0, frontal_area = 100.0)
        b = Building(; east = 80.0, north = 0.0, height = 30.0, frontal_area = 900.0)
        e = BuildingEnvelope([a, b])
        # Equal heights: the nearer governs, as NSR-23 §5 takes the nearest
        # building in the vicinity; the order of the list does not matter.
        @test equivalent_height(e) == 30.0
        @test equivalent_area(e) == 100.0
        @test equivalent_area(BuildingEnvelope([b, a])) == 100.0
        # A shed beside the stack does not cancel the wake of a reactor
        # building further out, which the 2021 inverse-distance mean did
        # (an equivalent 8.2 m).
        shed = Building(; east = 5.0, north = 0.0, height = 3.0, frontal_area = 10.0)
        reactor = Building(; east = 50.0, north = 0.0, height = 60.0, frontal_area = 3000.0)
        both = BuildingEnvelope([shed, reactor])
        @test equivalent_height(both) == 60.0
        @test equivalent_area(both) == 3000.0
        @test release_height(Site(; source = REFERENCE_STACK, atmosphere = REFERENCE_AIR,
            buildings = both,)) == 0
        # Outside three of its own heights a building plays no part.
        far = Building(; east = 200.0, north = 0.0, height = 60.0, frontal_area = 3000.0)
        @test equivalent_height(BuildingEnvelope([shed, far])) == 3.0
        # NSR-23 §5: in the vicinity when "smaller than 3 times" its height.
        edge = Building(; east = 180.0, north = 0.0, height = 60.0, frontal_area = 3000.0)
        inside = Building(; east = 179.0, north = 0.0, height = 60.0, frontal_area = 3000.0)
        @test equivalent_height(BuildingEnvelope([edge])) == 0
        @test equivalent_height(BuildingEnvelope([inside])) == 60.0
        # Equally tall and equally far: the larger frontal area, in either order.
        p = Building(; east = 20.0, north = 0.0, height = 30.0, frontal_area = 100.0)
        q = Building(; east = -20.0, north = 0.0, height = 30.0, frontal_area = 900.0)
        @test equivalent_area(BuildingEnvelope([p, q])) == 900.0
        @test equivalent_area(BuildingEnvelope([q, p])) == 900.0
        @test_throws ArgumentError BuildingEnvelope([p]; wake_coefficient = Inf)
        # A stack on the building, as in SRS-19 Example IV-3.
        roof = Building(; east = 0.0, north = 0.0, height = 30.0, frontal_area = 5000.0)
        @test equivalent_height(BuildingEnvelope([roof])) == 30.0
        @test_throws ArgumentError Building(;
            east = NaN, north = 0.0, height = 5.0, frontal_area = 1.0,)
        @test_throws ArgumentError Building(;
            east = 0.0, north = 0.0, height = Inf, frontal_area = 1.0,)
    end

    @testset "wake broadening by release height" begin
        e = BuildingEnvelope([
            Building(; east = 20.0, north = 0.0, height = 20.0, frontal_area = 400.0),
        ])
        σ = 10.0
        full = sqrt(σ^2 + DEFAULT_WAKE_COEFFICIENT * equivalent_area(e) / π)
        @test wake_broadened(σ, 0.0, e) ≈ full          # trapped in the cavity
        @test wake_broadened(σ, 19.0, e) ≈ full         # still below the tops
        @test wake_broadened(σ, 50.1, e) ≈ σ            # well clear
        between = wake_broadened(σ, 35.0, e)
        @test σ < between < full                        # interpolated
        # Setting the coefficient to zero disables the correction.
        off = BuildingEnvelope(
            [Building(; east = 20.0, north = 0.0, height = 20.0, frontal_area = 400.0)];
            wake_coefficient = 0.0,
        )
        @test wake_broadened(σ, 0.0, off) ≈ σ
        @test_throws DomainError wake_broadened(-1.0, 0.0, e)
        @test_throws DomainError wake_broadened(σ, -1.0, e)
    end

    # XOQDOQ Eqs. (8) and (11) and NSR-23 §5 limit the wake dilution of a
    # release at ground level to 3, which √3 on each parameter keeps. For
    # σ_y = 100 m, σ_z = 10 m and A = 3 × 10⁴ m² the vertical parameter binds
    # and the lateral does not: the centreline falls by 2.42, within the 3 of
    # the point form.
    @testset "the broadening limit of a ground-level release" begin
        big = BuildingEnvelope([
            Building(; east = 20.0, north = 0.0, height = 20.0, frontal_area = 3.0e4),
        ])
        Σy, Σz = wake_broadened(100.0, 0.0, big), wake_broadened(10.0, 0.0, big)
        @test Σz ≈ sqrt(3) * 10
        @test Σy ≈ sqrt(100.0^2 + 3.0e4 / π)
        @test Σy * Σz / (100 * 10) ≈ 2.4217 rtol = 1e-4
        unlimited = BuildingEnvelope(
            [Building(; east = 20.0, north = 0.0, height = 20.0, frontal_area = 3.0e4)];
            broadening_limit = Inf,)
        @test wake_broadened(10.0, 0.0, unlimited) ≈ sqrt(10.0^2 + 3.0e4 / π)
        # At σ = 0 the limit leaves zero and no limit the width of the wake
        # itself; neither is NaN.
        @test wake_broadened(0.0, 0.0, unlimited) ≈ sqrt(3.0e4 / π)
        @test wake_broadened(0.0, 0.0, big) == 0
        # An elevated plume is broadened without limit: a wider Σ_z brings it
        # down sooner, so the limit would lower its ground-level value. Below
        # the tops the full Eq. (9) applies, between the tops and 2.5 H_B the
        # interpolation of Eq. (10) on the unlimited value.
        for H in (1.0, 19.0)
            @test wake_broadened(10.0, H, big) ≈ sqrt(10.0^2 + 3.0e4 / π)
        end
        full = sqrt(10.0^2 + 3.0e4 / π)
        @test wake_broadened(10.0, 35.0, big) ≈ full - (35.0 - 20.0) / 30.0 * (full - 10.0)
        @test @inferred(wake_broadened(10.0f0, 0.0, big)) isa Float64
        @test @inferred(wake_broadened(10.0f0, 60.0, big)) isa Float64
    end
end

@testset "Building zones after SRS-19" begin
    e = BuildingEnvelope([
        Building(; east = 20.0, north = 0.0, height = 20.0, frontal_area = 400.0),
    ])
    # 2.5 H_B = 50 m and 2.5 √A_B = 50 m.
    @test CAVITY_EXTENT_FACTOR == 2.5
    @test building_zone(60.0, 10.0, e) === DISPLACEMENT_ZONE
    @test building_zone(50.0, 100.0, e) === WAKE_ZONE       # at the height bound, still disturbed
    @test building_zone(50.0, 50.0, e) === CAVITY_ZONE      # at the distance bound, still the cavity
    @test building_zone(0.0, 51.0, e) === WAKE_ZONE
    @test building_zone(0.0, 0.0, e) === CAVITY_ZONE
    @test building_zone(0.0, 1.0, BuildingEnvelope()) === DISPLACEMENT_ZONE
    @test_throws DomainError building_zone(-1.0, 10.0, e)
    @test_throws DomainError building_zone(10.0, -1.0, e)
end
