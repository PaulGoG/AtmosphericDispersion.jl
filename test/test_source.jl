@testset "Source and ambient state" begin
    # The stack and atmosphere of the 2021 run.
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

    @testset "fluxes" begin
        F = buoyancy_flux(stack, air)
        Fₘ = momentum_flux(stack, air)
        @test F ≈ (1.2 - 0.6) / 1.2 * STANDARD_GRAVITY * 10.0 * (2.33 / 2)^2
        @test Fₘ ≈ 0.6 / 1.2 * 10.0^2 * (2.33 / 2)^2
        @test F > 0        # lighter than air, so buoyant
        @test Fₘ > 0
        # A plume denser than the ambient air has a negative buoyancy flux.
        dense = StackSource(;
            height = 50.3,
            diameter = 2.33,
            exit_velocity = 10.0,
            exit_density = 2.0,
            exit_temperature = 324.0,
        )
        @test buoyancy_flux(dense, air) < 0
        @test_throws DomainError final_buoyant_rise(buoyancy_flux(dense, air), 5.0, 1e-3)
    end

    @testset "stability parameter" begin
        @test stability_parameter(air) > 0                     # inversion
        # Zero exactly at the dry adiabatic lapse rate.
        adiabatic = Atmosphere(;
            reference_speed = 4.0,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = -STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        @test stability_parameter(adiabatic) ≈ 0 atol = 1e-12
        # Negative in superadiabatic, genuinely unstable air.
        unstable = Atmosphere(;
            reference_speed = 4.0,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = -0.02,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        @test stability_parameter(unstable) < 0
        # The dry adiabatic lapse rate is about 9.8 K per km.
        @test STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT ≈ 0.00976 atol = 1e-5
    end

    @testset "validation" begin
        @test_throws ArgumentError StackSource(;
            height = 0,
            diameter = 2.33,
            exit_velocity = 10.0,
            exit_density = 0.6,
            exit_temperature = 324.0,
        )
        @test_throws ArgumentError StackSource(;
            height = 50.3,
            diameter = -1,
            exit_velocity = 10.0,
            exit_density = 0.6,
            exit_temperature = 324.0,
        )
        @test_throws ArgumentError StackSource(;
            height = 50.3,
            diameter = 2.33,
            exit_velocity = -1,
            exit_density = 0.6,
            exit_temperature = 324.0,
        )
        @test_throws ArgumentError Atmosphere(;
            reference_speed = -1,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = 0.0,
            surface = SURFACE_WATER,
            roughness = ROUGHNESS_PASTURE,
        )
        @test_throws ArgumentError Atmosphere(;
            reference_speed = 4.0,
            temperature = 0,
            density = 1.2,
            lapse_rate = 0.0,
            surface = SURFACE_WATER,
            roughness = ROUGHNESS_PASTURE,
        )
        @test_throws ArgumentError Atmosphere(;
            reference_speed = 4.0,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = NaN,
            surface = SURFACE_WATER,
            roughness = ROUGHNESS_PASTURE,
        )
    end
end
