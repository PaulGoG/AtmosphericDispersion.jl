@testset "Tabulated coefficients" begin
    @testset "lateral coefficient falls with stability" begin
        c = [lateral_coefficient(k) for k in PASQUILL_CLASSES]
        @test c == [0.22, 0.16, 0.11, 0.08, 0.06, 0.04]
        @test issorted(c; rev = true)
    end

    @testset "roughness lengths increase through the classes" begin
        z = [roughness_length(r) for r in ROUGHNESS_CLASSES]
        @test z == [0.01, 0.04, 0.1, 0.4, 1.0, 4.0]
        @test issorted(z)
    end

    @testset "profile exponent grows with roughness and stability" begin
        for s in WIND_PROFILE_SURFACES
            m = [profile_exponent(s, k) for k in PASQUILL_CLASSES]
            @test issorted(m)              # stability
            @test all(>(0), m)
        end
        for k in PASQUILL_CLASSES
            m = [profile_exponent(s, k) for s in WIND_PROFILE_SURFACES]
            @test issorted(m)              # roughness
        end
        @test profile_exponent(SURFACE_WATER, PASQUILL_A) == 0.03
        @test profile_exponent(SURFACE_FOREST_URBAN, PASQUILL_F) == 0.64
    end

    @testset "washout" begin
        for p in PRECIPITATION_TYPES
            Λ = [washout_coefficients(p, r) for r in PRECIPITATION_RATES]
            @test issorted([w.low for w in Λ])
            @test issorted([w.high for w in Λ])
            @test all(w -> w.low < w.high, Λ)
        end
        # Rain scavenges orders of magnitude faster than snow.
        for r in PRECIPITATION_RATES
            @test washout_coefficients(PRECIPITATION_RAIN, r).high >
                  100 * washout_coefficients(PRECIPITATION_SNOW, r).high
        end
        @test washout_coefficients(PRECIPITATION_RAIN, 1).high == 2.0e-4
        @test_throws ArgumentError washout_coefficients(PRECIPITATION_RAIN, 2.0)
    end
end
