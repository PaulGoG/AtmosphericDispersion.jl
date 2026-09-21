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

    @testset "contributions are weighted by inverse distance" begin
        a = Building(; east = 10.0, north = 0.0, height = 30.0, frontal_area = 100.0)
        b = Building(; east = 80.0, north = 0.0, height = 30.0, frontal_area = 900.0)
        e = BuildingEnvelope([a, b])
        # Equal heights, so the weighted height is that height regardless.
        @test equivalent_height(e) ≈ 30.0
        # The nearer, smaller building dominates the area.
        @test equivalent_area(e) ≈ (100 / 10 + 900 / 80) / (1 / 10 + 1 / 80)
        @test equivalent_area(e) < 900
        @test_throws ArgumentError BuildingEnvelope([
            Building(; east = 0.0, north = 0.0, height = 5.0, frontal_area = 1.0),
        ])
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
    end
end
