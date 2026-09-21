@testset "Plume rise" begin
    stack, air = REFERENCE_STACK, REFERENCE_AIR
    F = buoyancy_flux(stack, air)
    Fₘ = momentum_flux(stack, air)
    S = stability_parameter(air)

    @testset "transition distance" begin
        @test buoyancy_transition_distance(10.0) ≈ 14 * 10.0^(5 / 8)
        @test buoyancy_transition_distance(100.0) ≈ 34 * 100.0^(2 / 5)
        @test buoyancy_transition_distance(0.0) == 0
        # Continuous enough across the breakpoint to be a sane correlation.
        below = buoyancy_transition_distance(BUOYANCY_FLUX_BREAKPOINT - 1e-9)
        above = buoyancy_transition_distance(BUOYANCY_FLUX_BREAKPOINT)
        @test abs(below - above) / below < 0.1
        @test_throws DomainError buoyancy_transition_distance(-1.0)
    end

    @testset "rise grows with distance and saturates" begin
        for u in (1.0, 4.0, 10.0)
            xs = [1.0, 10.0, 100.0, 1000.0, 10_000.0]
            rises = [plume_rise(x, stack, air, u) for x in xs]
            @test issorted(rises)
            @test all(≥(0), rises)
            # Saturated well before the far field.
            @test plume_rise(1e5, stack, air, u) ≈ plume_rise(1e6, stack, air, u)
        end
    end

    @testset "a stronger wind bends the plume over" begin
        rises = [plume_rise(1000.0, stack, air, u) for u in (1.0, 2.0, 5.0, 10.0, 20.0)]
        @test issorted(rises; rev = true)
        @test all(>(0), rises)
    end

    @testset "stratification caps the rise" begin
        neutral = reference_atmosphere(;
            lapse_rate = -STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT)
        @test final_buoyant_rise(F, 4.0, stability_parameter(neutral)) ≥
              final_buoyant_rise(F, 4.0, S)
        # Unstable air has no stable ceiling, and must not raise a domain
        # error from a fractional power of a negative stability parameter.
        unstable = reference_atmosphere(; lapse_rate = -0.02)
        Su = stability_parameter(unstable)
        @test Su < 0
        @test isfinite(final_buoyant_rise(F, 4.0, Su))
        @test isfinite(final_momentum_rise(Fₘ, 10.0, 2.33, 4.0, Su))
        @test isfinite(plume_rise(1000.0, stack, unstable, 4.0))
    end

    @testset "domains" begin
        @test_throws DomainError final_buoyant_rise(F, 0.0, S)
        @test_throws DomainError final_momentum_rise(Fₘ, 10.0, 2.33, 0.0, S)
        @test_throws DomainError buoyant_rise(-1.0, F, 4.0, S)
        @test_throws DomainError momentum_rise(-1.0, Fₘ, 10.0, 2.33, 4.0, S)
    end

    # As for the dispersion parameters, reproduce the 2021 correlations
    # literally and require exact agreement wherever they are applicable,
    # i.e. in stably stratified air, which is the only case the original
    # could evaluate at all.
    @testset "fidelity to the 2021 implementation" begin
        original_X_0 = F -> F < 55 ? 14 * F^(5 / 8) : 34 * F^(2 / 5)
        original_hb_final = function (F, u, S)
            x_0 = original_X_0(F)
            min(
                2.6 * (F / (u * S))^(1 / 3),
                1.6 * F^(1 / 3) * (3.5 * x_0)^(2 / 3) / u,
                5.0 * F^(1 / 4) * S^(-3 / 8),
            )
        end
        original_hb = function (x, F, u, S)
            hbfinal = original_hb_final(F, u, S)
            hbtranzitie = 1.6 * F^(1 / 3) * x^(2 / 3) / u
            (x < 3.5 * original_X_0(F) && hbtranzitie <= hbfinal) ? hbtranzitie : hbfinal
        end
        original_hm_final = function (Fm, w_0, D, u, S)
            min(1.5 * w_0 * D / u, 4 * (Fm / S)^(1 / 4), 1.5 * (Fm / u)^(1 / 3) *
                                                         S^(-1 / 6))
        end
        original_hm = function (x, Fm, w_0, D, u, S)
            hmfinal = original_hm_final(Fm, w_0, D, u, S)
            hmtranzitie = 1.89 * (w_0^2 * D / (u * (w_0 + 3u)))^(2 / 3) * x^(1 / 3)
            hmtranzitie <= hmfinal ? hmtranzitie : hmfinal
        end
        original_hmb = function (x, F, Fm, w_0, D, u, S)
            hmbfinal = original_hm_final(Fm, w_0, D, u, S) + original_hb_final(F, u, S)
            hmbtranzitie = 3^(1 / 3) *
                           (Fm * x / ((1 / 3 + u / w_0)^2 * u^2) + F * x^2 / (0.5 * u^3))^(1 /
                                                                                           3)
            hmbtranzitie <= hmbfinal ? hmbtranzitie : hmbfinal
        end
        original_rise = function (x, F, Fm, w_0, D, u, S)
            hbfinal = original_hb_final(F, u, S)
            hmfinal = original_hm_final(Fm, w_0, D, u, S)
            if abs(hbfinal - hmfinal) * 2 / (hbfinal + hmfinal) <= 0.1
                return original_hmb(x, F, Fm, w_0, D, u, S)
            elseif hmfinal > hbfinal
                return original_hm(x, Fm, w_0, D, u, S)
            else
                return original_hb(x, F, u, S)
            end
        end

        w₀, D = 10.0, 2.33
        for u in (0.5, 1.0, 4.0, 7.3, 15.0), Sv in (1e-4, 8.63e-4, 1e-3, 5e-3)

            @test final_buoyant_rise(F, u, Sv, NSR23_RISE) == original_hb_final(F, u, Sv)
            @test final_momentum_rise(Fₘ, w₀, D, u, Sv, NSR23_RISE) ==
                  original_hm_final(Fₘ, w₀, D, u, Sv)
            for x in (1.0, 10.0, 137.0, 1000.0, 10_000.0)
                @test buoyant_rise(x, F, u, Sv, NSR23_RISE) == original_hb(x, F, u, Sv)
                @test momentum_rise(x, Fₘ, w₀, D, u, Sv, NSR23_RISE) ==
                      original_hm(x, Fₘ, w₀, D, u, Sv)
                @test combined_rise(x, F, Fₘ, w₀, u, Sv, D, NSR23_RISE) ==
                      original_hmb(x, F, Fₘ, w₀, D, u, Sv)
            end
        end
    end

    @testset "both methods of plume_rise honour the coefficients" begin
        # A cold jet in a neutral atmosphere: no buoyancy, no stratification,
        # so the final rise is c w₀D/u and the two sets differ by exactly 3/1.5.
        jet = StackSource(;
            height = 30.0,
            diameter = 1.0,
            exit_velocity = 20.0,
            exit_density = 1.2,
            exit_temperature = 287.0,
        )
        neutral = reference_atmosphere(; lapse_rate = -0.0098)
        @test buoyancy_flux(jet, neutral) == 0
        u, x = 5.0, 1e5
        @test plume_rise(x, jet, neutral, u) == 3 * 20.0 * 1.0 / u
        @test plume_rise(x, jet, neutral, u, NSR23_RISE) == 1.5 * 20.0 * 1.0 / u
        for rise in (BRIGGS_RISE, XOQDOQ_RISE, NSR23_RISE)
            site = Site(; source = jet, atmosphere = neutral, rise)
            v = transport_wind_speed(site, PASQUILL_D)
            @test plume_rise(x, site, PASQUILL_D) == plume_rise(x, jet, neutral, v, rise)
        end
    end
end
