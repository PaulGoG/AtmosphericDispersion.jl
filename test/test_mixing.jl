# HPA-RPD-058 §3.2.2.1 Eqs. (3.4) and (3.5), and Table 3.5(a), which
# the report attributes to Clarke (1979) and Jones (1980).
@testset "the mixing layer is HPA-RPD-058" begin
    published = (1300.0, 900.0, 850.0, 800.0, 400.0, 100.0)
    for (i, class) in enumerate(PASQUILL_CLASSES)
        @test mixing_depth(class) == published[i]
        @test mixing_depth(class, MIXING_UNBOUNDED) == Inf
        @test mixing_depth(class, MIXING_UNIFORM_800) == 800.0
    end
    # deeper in unstable air, shallower in stable, monotonic past A
    @test issorted(published[2:end]; rev = true)
    @test RECOMMENDED_MIXING_DEPTH == 800.0

    # Without a lid the factor is the ordinary ground-reflected Gaussian.
    for H in (20.0, 100.0), Σz in (10.0, 60.0, 400.0), z in (0.0, 50.0)
        plain = (exp(-(z - H)^2 / (2Σz^2)) + exp(-(z + H)^2 / (2Σz^2))) / (sqrt(2π) * Σz)
        @test vertical_factor(z, H, Σz, Inf) ≈ plain
    end

    # Eq. (3.5): once σ_z reaches the depth the profile is uniform.
    for A in (400.0, 800.0, 1300.0), H in (10.0, 100.0)
        @test vertical_factor(0.0, H, A, A) ≈ 1 / A
        @test vertical_factor(A / 2, H, 3A, A) ≈ 1 / A
    end

    # Nothing crosses the lid.
    for A in (100.0, 800.0)
        @test vertical_factor(A + 1, 50.0, 30.0, A) == 0
        @test vertical_factor(A + 1, 50.0, 3A, A) == 0
    end

    # The released activity is conserved under the lid, which is the
    # statement the image sum has to satisfy and the truncation at
    # |s| = 1 could have broken. All of these have H < A, which is the
    # regime the lid model covers.
    for (H, Σz, A) in (
        (100.0, 30.0, 800.0),
        (100.0, 300.0, 800.0),
        (50.0, 900.0, 800.0),
        (20.0, 40.0, 100.0),
        (10.0, 250.0, 100.0),
    )
        mass, _ = quadgk(z -> vertical_factor(z, H, Σz, A), 0.0, A; rtol = 1e-12)
        @test mass ≈ 1 rtol = 1e-8
    end
    # Those cases sit below Σz = 0.4 A or in the uniform branch, where
    # conservation is exact. Between them the truncation at |s| = 1
    # loses activity, most just below the switch; it is bounded here so
    # that it cannot grow unnoticed.
    A = 800.0
    conserved(H, ratio) = quadgk(z -> vertical_factor(z, H, ratio * A, A), 0.0, A)[1]
    for H = 0.0:100.0:A
        @test 1 - 5e-4 < conserved(H, 0.6) ≤ 1 + 1e-12
        @test 0.977 < conserved(H, 0.999) ≤ 1 + 1e-12
    end
    @test conserved(80.0, 0.999) ≈ 0.997 atol = 1e-3
    @test conserved(A, 0.999) ≈ 0.977 atol = 1e-3
    # And the step across the switch is what the docstring states.
    step(H, z) = abs(vertical_factor(z, H, prevfloat(A), A) * A - 1)
    @test step(0.0, 0.0) ≈ 0.014 atol = 1e-3
    @test maximum(step(H, z) for H = 0.0:100.0:A, z = 0.0:100.0:A) ≈ 0.040 atol = 1e-3

    # A lid can only raise the ground-level concentration: it reflects
    # material back down that would otherwise have gone on rising.
    for H in (20.0, 100.0), Σz in (50.0, 200.0, 700.0), A in (400.0, 800.0)
        @test vertical_factor(0.0, H, Σz, A) ≥ vertical_factor(0.0, H, Σz, Inf)
    end

    # A release strictly above the lid is not trapped by it. HPA's
    # Diagram 3.1 places the source below the inversion; a plume above
    # one is decoupled until the inversion breaks, which is a different
    # model. The reference case meets this in class F: 100 m depth
    # against a 103 m release. A release exactly at the lid is still
    # capped — Table 3.7 tabulates one and only the capped form
    # reproduces it.
    @test vertical_factor(0.0, 100.0, 50.0, 100.0) != vertical_factor(0.0, 100.0, 50.0, Inf)
    for H in (150.0, 101.0), Σz in (20.0, 200.0)
        @test vertical_factor(0.0, H, Σz, 100.0) == vertical_factor(0.0, H, Σz, Inf)
        @test vertical_factor(500.0, H, Σz, 100.0) == vertical_factor(500.0, H, Σz, Inf)
    end
end
