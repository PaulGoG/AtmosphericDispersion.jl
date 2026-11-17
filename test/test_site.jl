@testset "Site" begin
    stack, air, site = REFERENCE_STACK, REFERENCE_AIR, REFERENCE_SITE
    tall = BuildingEnvelope([
        Building(; east = 50.0, north = 0.0, height = 100.0, frontal_area = 2000.0),
    ])

    @testset "precomputation matches the direct evaluation" begin
        @test release_height(site) == wake_height(stack, air, BuildingEnvelope())
        @test site.buoyancy == buoyancy_flux(stack, air)
        @test site.momentum == momentum_flux(stack, air)
        @test site.stability == stability_parameter(air)
    end

    @testset "downwash" begin
        # A vigorous efflux in a light wind is not drawn down.
        @test downwash_height(stack, air) == stack.height
        # A weak efflux in a strong wind is drawn down.
        calm = reference_stack(; exit_velocity = 1.0)
        windy = reference_atmosphere(; reference_speed = 15.0)
        @test downwash_height(calm, windy) < calm.height
        # In calm air the test takes the 1 m/s floor: 2 (1.5 − 1/1) D lower.
        still = reference_atmosphere(; reference_speed = 0.0)
        @test downwash_height(calm, still) ≈ 50.3 - 2 * 0.5 * 2.33
        # A short, wide stack with no efflux is drawn down to the ground and
        # no further, and with no building it is released there.
        for h in (3.0, 2.0)
            stub = reference_stack(; height = h, diameter = 1.0, exit_velocity = 0.0)
            @test downwash_height(stub, still) == 0
            @test wake_height(stub, still, BuildingEnvelope()) == 0
        end
        # A trickle of efflux gives a stable rise of micrometres; the mean wind
        # over so thin a layer is the wind at its base.
        trickle = Site(; source = reference_stack(; exit_velocity = 1e-20), atmosphere = air)
        for class in STABLE_CLASSES
            @test stable_rise_wind_speed(trickle, class) ≈
                  max(
                wind_speed(air.reference_speed, release_height(trickle), air.surface,
                    class,),
                MINIMUM_WIND_SPEED,) rtol = 1e-6
        end
    end

    @testset "wake height" begin
        # No buildings: release height is the downwash-corrected stack height.
        @test wake_height(stack, air, BuildingEnvelope()) == downwash_height(stack, air)
        # A building taller than the stack traps the plume at ground level.
        @test wake_height(stack, air, tall) == 0
        # Between one and two and a half building heights: entrained, by
        # 1.5 h − 0.6 H₁ (NSR-23 Eq. (14)), only when the mean wind at the stack
        # reaches 5 m/s.
        mid = BuildingEnvelope([
            Building(; east = 30.0, north = 0.0, height = 30.0, frontal_area = 900.0),
        ])
        light = reference_atmosphere(; reference_speed = 2.0)
        windy = reference_atmosphere(; reference_speed = 8.0)
        for a in (light, windy)
            H₁ = downwash_height(stack, a)
            @test 30 ≤ H₁ < 75
            expected = mean_wind_speed(a, H₁) < WAKE_WIND_THRESHOLD ? H₁ :
                       H₁ - (1.5 * 30 - 0.6 * H₁)
            @test wake_height(stack, a, mid) ≈ expected
        end
        @test mean_wind_speed(light, downwash_height(stack, light)) < WAKE_WIND_THRESHOLD
        @test mean_wind_speed(windy, downwash_height(stack, windy)) ≥ WAKE_WIND_THRESHOLD
        # At 2.5 h the correction vanishes, so the boundary is continuous.
        H₁ = downwash_height(stack, windy)
        edge = BuildingEnvelope([
            Building(; east = 30.0, north = 0.0, height = H₁ / 2.5, frontal_area = 900.0),
        ])
        @test wake_height(stack, windy, edge) ≈ H₁
    end

    @testset "transport wind is floored at the reference height" begin
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

    @testset "the mixing layer of a site" begin
        @test mixing_layer(site) === MIXING_TABULATED        # the default
        @test mixing_layer(UNBOUNDED_SITE) === MIXING_UNBOUNDED
    end
end

@testset "PrescribedPlume" begin
    site = REFERENCE_SITE
    distances = (10.0, 1000.0, 1e5)

    @testset "validation" begin
        @test_throws ArgumentError PrescribedPlume(site; height = -1.0)
        for v in (0.0, -1.0)
            @test_throws ArgumentError PrescribedPlume(site; wind = v)
            @test_throws ArgumentError PrescribedPlume(site; lateral = v)
            @test_throws ArgumentError PrescribedPlume(site; vertical = v)
        end
        # A ground-level release is a legitimate prescription.
        @test effective_height(100.0, PrescribedPlume(site; height = 0), PASQUILL_D) == 0
    end

    @testset "what is not prescribed comes from the site" begin
        bare = PrescribedPlume(site)
        @test mixing_layer(bare) === mixing_layer(site)
        for class in PASQUILL_CLASSES
            @test transport_wind_speed(bare, class) == transport_wind_speed(site, class)
            for x in distances
                @test effective_height(x, bare, class) == effective_height(x, site, class)
                @test corrected_lateral_dispersion(x, bare, class) ==
                      corrected_lateral_dispersion(x, site, class)
                @test corrected_lateral_dispersion(
                    x,
                    bare,
                    class;
                    release_duration = 3600,
                ) == corrected_lateral_dispersion(x, site, class; release_duration = 3600)
                @test corrected_vertical_dispersion(x, bare, class) ==
                      corrected_vertical_dispersion(x, site, class)
            end
        end
    end

    @testset "what is prescribed is returned as given" begin
        H, u, σy, σz = 150.0, 4.0, 157.0, 110.0
        plume = PrescribedPlume(site; height = H, wind = u, lateral = σy, vertical = σz)
        @test mixing_layer(plume) === mixing_layer(site)
        for class in PASQUILL_CLASSES
            @test transport_wind_speed(plume, class) == u
            for x in distances
                @test effective_height(x, plume, class) == H
                @test corrected_lateral_dispersion(x, plume, class) == σy
                @test corrected_vertical_dispersion(x, plume, class) == σz
            end
        end
        # One field at a time: the others are untouched. There are no buildings
        # here, so a prescribed height does not affect the dispersion parameters.
        only_wind = PrescribedPlume(site; wind = u)
        only_height = PrescribedPlume(site; height = H)
        for class in PASQUILL_CLASSES, x in distances

            @test effective_height(x, only_wind, class) == effective_height(x, site, class)
            @test transport_wind_speed(only_height, class) ==
                  transport_wind_speed(site, class)
            @test corrected_vertical_dispersion(x, only_height, class) ==
                  corrected_vertical_dispersion(x, site, class)
        end
    end

    # A constant σ_z is a statement about one distance, and the integral runs
    # over all of them.
    @testset "the depletion integral refuses a prescribed σ_z" begin
        @test_throws ArgumentError depletion_integral(
            1000.0,
            PrescribedPlume(site; vertical = 110.0),
            PASQUILL_D,
        )
        @test_throws ArgumentError dry_depletion_factor(
            1000.0,
            PrescribedPlume(site; vertical = 110.0),
            PASQUILL_D,
            TRITIATED_WATER,
        )
        @test depletion_integral(
            1000.0,
            PrescribedPlume(site; height = 80.0, wind = 5.0, lateral = 157.0),
            PASQUILL_D,
        ) > 0
    end

    # The wake broadening is decided by the height of the plume, and a plume
    # whose height is prescribed is judged on that height.
    @testset "the wake test sees the prescribed height" begin
        envelope = BuildingEnvelope([
            Building(; east = 30.0, north = 0.0, height = 40.0, frontal_area = 1600.0),
        ])
        built = Site(;
            source = REFERENCE_STACK,
            atmosphere = REFERENCE_AIR,
            buildings = envelope,
        )
        h = equivalent_height(envelope)
        clear = PrescribedPlume(built; height = 2.5h + 1)
        inside = PrescribedPlume(built; height = h / 2)
        roughness = REFERENCE_AIR.roughness
        for class in PASQUILL_CLASSES, x in (100.0, 500.0)

            σy = lateral_dispersion(x, class)
            σz = vertical_dispersion(x, class, roughness)
            # The site's own plume is still within reach of the wake in the
            # stable classes, whose rise is the lowest.
            if class in STABLE_CLASSES
                @test effective_height(x, built, class) < 2.5h
                @test corrected_vertical_dispersion(x, built, class) > σz
            end
            @test corrected_vertical_dispersion(x, clear, class) == σz
            @test corrected_lateral_dispersion(x, clear, class) == σy
            @test corrected_vertical_dispersion(x, inside, class) ==
                  wake_broadened(σz, h / 2, envelope)
            @test corrected_lateral_dispersion(x, inside, class) ==
                  wake_broadened(σy, h / 2, envelope)
        end
    end
end

@testset "The wind of the stable rise of a site" begin
    site = REFERENCE_SITE       # stably stratified and buoyant: the mean applies
    at_release = Site(;
        source = REFERENCE_STACK,
        atmosphere = REFERENCE_AIR,
        stable_rise_wind = WIND_AT_RELEASE_HEIGHT,
    )
    @test site.stable_rise_wind === WIND_MEAN_OVER_RISE
    c, F, S = BRIGGS_RISE.stable_final, site.buoyancy, site.stability
    z₀ = max(release_height(site), REFERENCE_HEIGHT)
    # The stable limits apply in E and F only: in A to D the mean is not
    # formed and the plume takes the neutral limit.
    for class in PASQUILL_CLASSES
        class in STABLE_CLASSES && continue
        @test stability_parameter(site, class) == 0
        @test stable_rise_wind_speed(site, class) == transport_wind_speed(site, class)
        u = transport_wind_speed(site, class)
        x₀ = buoyancy_transition_distance(F)
        @test plume_rise(1e5, site, class) ≈ 1.6 * F^(1 / 3) * (3.5 * x₀)^(2 / 3) / u
    end
    for class in STABLE_CLASSES
        @test stability_parameter(site, class) == S
        u = transport_wind_speed(site, class)
        ū = stable_rise_wind_speed(site, class)
        @test stable_rise_wind_speed(at_release, class) == u
        # Positive shear: the mean over the rise exceeds the wind at its base.
        @test ū > u
        # It is the mean of the profile over exactly the rise it produces.
        Δh = c * (F / (ū * S))^(1 / 3)
        @test layer_mean_wind_speed(
            REFERENCE_AIR.reference_speed, z₀, z₀ + Δh, REFERENCE_AIR.surface, class,) ≈ ū rtol = 1e-10
        # A faster wind cannot raise the plume.
        @test effective_height(1e5, site, class) ≤ effective_height(1e5, at_release, class)
    end
    # In class F of the reference case the stable limit binds, and the mean
    # wind lowers it.
    @test effective_height(1e5, site, PASQUILL_F) <
          effective_height(1e5, at_release, PASQUILL_F)

    # No stable rise to average over: unstable air, or a plume without buoyancy.
    unstable = Site(;
        source = REFERENCE_STACK,
        atmosphere = reference_atmosphere(; lapse_rate = -0.02),
    )
    jet = Site(; source = reference_stack(; exit_density = 1.2), atmosphere = REFERENCE_AIR)
    @test jet.buoyancy == 0
    for class in PASQUILL_CLASSES
        @test stable_rise_wind_speed(unstable, class) ==
              transport_wind_speed(unstable, class)
        @test stable_rise_wind_speed(jet, class) == transport_wind_speed(jet, class)
    end
    # Calm air: every wind of the site is held at 1 m/s (EPA ISC3 §1.1.3),
    # where the neutral rise, which varies as 1/u, would otherwise be
    # unbounded.
    calm = Site(;
        source = REFERENCE_STACK,
        atmosphere = reference_atmosphere(; reference_speed = 0.0),
    )
    for class in PASQUILL_CLASSES
        @test transport_wind_speed(calm, class) == MINIMUM_WIND_SPEED
        @test stable_rise_wind_speed(calm, class) == MINIMUM_WIND_SPEED
        @test isfinite(effective_height(1e5, calm, class))
    end
end

@testset "The dispersion scheme of a site" begin
    @test dispersion_scheme(REFERENCE_SITE) === DISPERSION_HOSKER
    for s in DISPERSION_SCHEMES
        site = Site(; source = REFERENCE_STACK, atmosphere = REFERENCE_AIR, dispersion = s)
        @test dispersion_scheme(site) === s
        @test dispersion_scheme(PrescribedPlume(site)) === s
        for class in PASQUILL_CLASSES, x in (100.0, 1000.0, 10_000.0)

            @test corrected_lateral_dispersion(x, site, class) ==
                  lateral_dispersion(x, class, s)
            @test corrected_lateral_dispersion(x, site, class; release_duration = 3600) ==
                  lateral_dispersion(x, class, s; release_duration = 3600)
            @test corrected_vertical_dispersion(x, site, class) ==
                  vertical_dispersion(x, class, s, REFERENCE_AIR.roughness)
            # With no buildings a prescribed height changes neither parameter.
            plume = PrescribedPlume(site; height = 10.0)
            @test corrected_vertical_dispersion(x, plume, class) ==
                  corrected_vertical_dispersion(x, site, class)
            @test corrected_lateral_dispersion(x, plume, class) ==
                  corrected_lateral_dispersion(x, site, class)
        end
    end
end
