@testset "Site" begin
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

    @testset "precomputation matches the direct evaluation" begin
        @test release_height(site) == wake_height(stack, air, BuildingEnvelope())
        @test site.buoyancy == buoyancy_flux(stack, air)
        @test site.momentum == momentum_flux(stack, air)
        @test site.stability == stability_parameter(air)
    end

    @testset "downwash" begin
        # A vigorous efflux in a light wind clears the stack untouched.
        @test downwash_height(stack, air) == stack.height
        # A weak efflux in a strong wind is drawn down.
        calm = StackSource(;
            height = 50.3,
            diameter = 2.33,
            exit_velocity = 1.0,
            exit_density = 0.6,
            exit_temperature = 324.0,
        )
        windy = Atmosphere(;
            reference_speed = 15.0,
            temperature = 287.0,
            density = 1.2,
            lapse_rate = 2e-2,
            surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        @test downwash_height(calm, windy) < calm.height
    end

    @testset "wake height" begin
        # No buildings: release height is the downwash-corrected stack height.
        @test wake_height(stack, air, BuildingEnvelope()) == downwash_height(stack, air)
        # A building taller than the stack traps the plume at ground level.
        tall = BuildingEnvelope([
            Building(; east = 50.0, north = 0.0, height = 100.0, frontal_area = 2000.0),
        ])
        @test wake_height(stack, air, tall) == 0
    end

    @testset "transport wind is floored at the reference height" begin
        tall = BuildingEnvelope([
            Building(; east = 50.0, north = 0.0, height = 100.0, frontal_area = 2000.0),
        ])
        trapped = Site(; source = stack, atmosphere = air, buildings = tall)
        @test release_height(trapped) == 0
        # Without the floor this would be zero and every dilution factor
        # would divide by it.
        for class in PASQUILL_CLASSES
            u = transport_wind_speed(trapped, class)
            @test u > 0
            @test u == wind_speed(4.0, REFERENCE_HEIGHT, SURFACE_AGRICULTURAL, class)
        end
    end

    @testset "effective height rises with distance" begin
        for class in PASQUILL_CLASSES
            hs = [effective_height(x, site, class) for x in (1.0, 100.0, 1000.0, 1e5)]
            @test issorted(hs)
            @test all(≥(release_height(site)), hs)
        end
    end
end
