@testset "Source and ambient state" begin
    stack, air = REFERENCE_STACK, REFERENCE_AIR

    @testset "fluxes" begin
        F = buoyancy_flux(stack, air)
        Fₘ = momentum_flux(stack, air)
        @test F ≈ (1.2 - 0.6) / 1.2 * STANDARD_GRAVITY * 10.0 * (2.33 / 2)^2
        @test Fₘ ≈ 0.6 / 1.2 * 10.0^2 * (2.33 / 2)^2
        @test F > 0        # lighter than air, so buoyant
        @test Fₘ > 0
        # A plume denser than the ambient air has a negative buoyancy flux.
        dense = reference_stack(; exit_density = 2.0)
        @test buoyancy_flux(dense, air) < 0
        @test_throws DomainError final_buoyant_rise(buoyancy_flux(dense, air), 5.0, 1e-3)
    end

    @testset "stability parameter" begin
        @test stability_parameter(air) > 0                     # inversion
        # Zero exactly at the dry adiabatic lapse rate.
        adiabatic = reference_atmosphere(;
            lapse_rate = -STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT)
        @test stability_parameter(adiabatic) ≈ 0 atol = 1e-12
        # Negative in superadiabatic, genuinely unstable air.
        unstable = reference_atmosphere(; lapse_rate = -0.02)
        @test stability_parameter(unstable) < 0
        # The dry adiabatic lapse rate is about 9.8 K per km.
        @test STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT ≈ 0.00976 atol = 1e-5
    end

    @testset "validation" begin
        @test_throws ArgumentError reference_stack(; height = 0)
        @test_throws ArgumentError reference_stack(; diameter = -1)
        @test_throws ArgumentError reference_stack(; exit_velocity = -1)
        @test_throws ArgumentError reference_atmosphere(; reference_speed = -1)
        @test_throws ArgumentError reference_atmosphere(; temperature = 0)
        @test_throws ArgumentError reference_atmosphere(; lapse_rate = NaN)
    end
end
