# HPA-RPD-058 §3.2.2.1 Eqs. (3.4) and (3.5), and Table 3.5(a), which the report
# attributes to Clarke (1979) and Jones (1980). The reference values here are
# computed from the published series, never from the package.
@testset "Mixing layer" begin
    # Eq. (3.4) summed by brute force far past convergence. At Σ_z = 5A, the
    # widest case below, the first neglected image is down by exp(−288).
    image_sum(z, H, Σz, A) = sum(
        exp(-(z - H - 2s * A)^2 / (2Σz^2)) + exp(-(z + H - 2s * A)^2 / (2Σz^2))
    for
    s in -60:60
    ) / (sqrt(2π) * Σz)
    # Σ_z/A on both sides of the change of representation at 0.7, and hard
    # against it.
    ratios = (0.05, 0.3, 0.69, 0.6999999, 0.7, 0.71, 1.0, 2.0, 5.0)

    @testset "MixingLayer" begin
        # Table 3.5(a), entry by entry.
        published = (1300.0, 900.0, 850.0, 800.0, 400.0, 100.0)
        @test MIXING_TABULATED.depths == published
        @test MIXING_TABULATED.above_lid === RISE_INHIBITED
        @test RECOMMENDED_MIXING_DEPTH == 800.0
        @test MixingLayer(800.0).depths == ntuple(_ -> 800.0, 6)
        recommended = MixingLayer(RECOMMENDED_MIXING_DEPTH)
        for (i, class) in enumerate(PASQUILL_CLASSES)
            @test mixing_depth(class) == published[i]
            @test mixing_depth(class, MIXING_UNBOUNDED) == Inf
            @test mixing_depth(class, recommended) == 800.0
        end
        @test MixingLayer(published; above_lid = FULL_PENETRATION).above_lid ===
              FULL_PENETRATION
        # Integer depths are as good as floating-point ones.
        @test MixingLayer((1300, 900, 850, 800, 400, 100)) == MIXING_TABULATED

        @test_throws ArgumentError MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0))
        @test_throws ArgumentError MixingLayer(fill(800.0, 7))
        @test_throws ArgumentError MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0, 0.0))
        @test_throws ArgumentError MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0, -100.0))
        @test_throws ArgumentError MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0, NaN))
        @test_throws ArgumentError MixingLayer((1300.0, 900.0, 850.0, 800.0, 400.0, "100"))
        @test_throws ArgumentError MixingLayer(0.0)
        @test_throws ArgumentError MixingLayer(-800.0)
    end

    @testset "domains" begin
        @test_throws DomainError vertical_factor(0.0, 50.0, 0.0, 800.0)
        @test_throws DomainError vertical_factor(0.0, 50.0, 30.0, 0.0)
        @test_throws DomainError vertical_factor(-1.0, 50.0, 30.0, 800.0)
        @test_throws DomainError vertical_factor(0.0, -1.0, 30.0, 800.0)
    end

    @testset "the profile is the image sum of Eq. (3.4)" begin
        A = 800.0
        for ratio in ratios, H in (0.0, A / 10, A / 2, A), z in (0.0, 123.0, A)
            Σz = ratio * A
            @test vertical_factor(z, H, Σz, A) ≈ image_sum(z, H, Σz, A) rtol = 1e-12
        end
    end

    # What the lid has to satisfy whatever the representation: the released
    # activity stays between the ground and the lid, a source at the lid
    # included.
    @testset "activity under the lid is conserved" begin
        for A in (800.0, 100.0), ratio in ratios, H in (0.0, A / 10, A / 2, A)
            Σz = ratio * A
            activity, _ = quadgk(z -> vertical_factor(z, H, Σz, A), 0, A; rtol = 1e-12)
            @test activity ≈ 1 atol = 1e-10
        end
    end

    @testset "limits" begin
        # The cosine series at Σ_z = A, source and receptor on the ground: every
        # harmonic enters with weight one, and the third is below round-off.
        for A in (100.0, 800.0)
            @test A * vertical_factor(0.0, 0.0, A, A) ≈ 1 + 2exp(-π^2 / 2) + 2exp(-2π^2) rtol = 1e-12
        end
        # Eq. (3.5): a plume that has filled the layer is uniform over it.
        for A in (100.0, 800.0), H in (0.0, A / 10, A / 2, A), z in (0.0, A / 3, A)
            @test A * vertical_factor(z, H, 3A, A) ≈ 1 atol = 1e-12
        end
        # Without a lid the factor is the ordinary ground-reflected Gaussian.
        for H in (20.0, 100.0), Σz in (10.0, 60.0, 400.0), z in (0.0, 50.0)
            plain = (exp(-(z - H)^2 / (2Σz^2)) + exp(-(z + H)^2 / (2Σz^2))) /
                    (sqrt(2π) * Σz)
            @test vertical_factor(z, H, Σz, Inf) == plain
            @test vertical_factor(z, H, Σz, Inf, FULL_PENETRATION) == plain
        end
        # Nothing crosses the lid.
        for A in (100.0, 800.0), Σz in (30.0, 3A)

            @test vertical_factor(A + 1, 50.0, Σz, A) == 0
            @test vertical_factor(nextfloat(A), 50.0, Σz, A) == 0
        end
    end

    @testset "the crosswind-integrated factor without a lid" begin
        for H in (20.0, 100.0), Σz in (10.0, 60.0, 400.0)

            @test crosswind_integrated_factor(H, Σz, Inf) ≈
                  sqrt(2 / π) * exp(-H^2 / (2Σz^2)) / Σz rtol = 1e-14
        end
    end

    # A lid can only raise the ground-level concentration: it reflects material
    # back down that would otherwise have gone on rising.
    @testset "a lid raises the ground-level factor" begin
        for H in (20.0, 100.0), Σz in (50.0, 200.0, 700.0), A in (400.0, 800.0)
            @test vertical_factor(0.0, H, Σz, A) ≥ vertical_factor(0.0, H, Σz, Inf)
        end
    end

    # NRPB-R157 §B2.3: the capping inversion arrests the rise, so a plume that
    # would have gone above the lid is dispersed from it. The reference case
    # meets this in class F, 100 m of depth against a release near 103 m.
    @testset "a plume above the lid, rise inhibited" begin
        A = 800.0
        for Σz in (80.0, 240.0, 800.0, 2400.0), z in (0.0, 123.0, A)

            at_lid = vertical_factor(z, A, Σz, A)
            @test at_lid > 0
            for H in (nextfloat(A), A + 0.1, 1.5A, 10A)
                @test vertical_factor(z, H, Σz, A) == at_lid
                @test vertical_factor(z, H, Σz, A, RISE_INHIBITED) == at_lid
            end
            # Continuous across H = A: the plume held at the lid is the limit of
            # the plume just below it.
            below = vertical_factor(z, A - 0.1, Σz, A)
            above = vertical_factor(z, A + 0.1, Σz, A)
            @test abs(above - below) < 1e-3 * below
        end
    end

    # EPA ISC3 User's Guide vol. II §1.1.6.1: a plume whose effective height
    # exceeds the mixing height leaves the layer and the ground sees nothing.
    @testset "a plume above the lid, full penetration" begin
        A = 800.0
        for Σz in (80.0, 240.0, 800.0, 2400.0), z in (0.0, 123.0, A)

            for H in (nextfloat(A), A + 0.1, 1.5A, 10A)
                @test vertical_factor(z, H, Σz, A, FULL_PENETRATION) == 0
            end
            # At and below the lid the rule is not consulted.
            for H in (0.0, A / 10, A / 2, prevfloat(A), A)
                @test vertical_factor(z, H, Σz, A, FULL_PENETRATION) ==
                      vertical_factor(z, H, Σz, A, RISE_INHIBITED)
            end
        end
    end

    # ISC3 applies its mixing heights to the unstable and neutral categories and
    # takes stable air as unbounded.
    @testset "an ISC3-style layer" begin
        isc3 = MixingLayer((1500, 1200, 1000, 900, Inf, Inf); above_lid = FULL_PENETRATION)
        for Σz in (50.0, 400.0), z in (0.0, 60.0)

            for class in (PASQUILL_E, PASQUILL_F), H in (100.0, 2000.0)

                @test vertical_factor(z, H, Σz, isc3, class) ==
                      vertical_factor(z, H, Σz, Inf)
                @test crosswind_integrated_factor(H, Σz, isc3, class) ==
                      vertical_factor(0.0, H, Σz, Inf)
            end
            for class in (PASQUILL_A, PASQUILL_B, PASQUILL_C, PASQUILL_D)
                A = mixing_depth(class, isc3)
                @test vertical_factor(z, A + 1, Σz, isc3, class) == 0
                @test crosswind_integrated_factor(A + 1, Σz, isc3, class) == 0
                # Below the lid the class method is the depth method.
                @test vertical_factor(z, 100.0, Σz, isc3, class) ==
                      vertical_factor(z, 100.0, Σz, A, FULL_PENETRATION)
            end
        end
        # The same depths with the rise inhibited hold the plume at the lid.
        held = MixingLayer(isc3.depths)
        @test vertical_factor(0.0, 2000.0, 400.0, held, PASQUILL_D) ==
              vertical_factor(0.0, 900.0, 400.0, 900.0)
    end
end
