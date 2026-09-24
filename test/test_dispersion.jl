@testset "Dispersion parameters" begin
    @testset "wind profile" begin
        u₁₀ = 4.0
        for s in WIND_PROFILE_SURFACES, k in PASQUILL_CLASSES

            @test wind_speed(u₁₀, 10, s, k) ≈ u₁₀       # the reference height
            @test wind_speed(u₁₀, 50, s, k) > u₁₀       # shear
            @test wind_speed(u₁₀, 5, s, k) < u₁₀
            # Held constant above the ceiling.
            @test wind_speed(u₁₀, 200, s, k) == wind_speed(u₁₀, 1000, s, k)
            @test wind_speed(u₁₀, 201, s, k) == wind_speed(u₁₀, 200, s, k)
            # Monotone below it.
            @test wind_speed(u₁₀, 20, s, k) < wind_speed(u₁₀, 100, s, k)
            @test wind_speed(0, 50, s, k) == 0
        end
        @test_throws DomainError wind_speed(4.0, 0, SURFACE_WATER, PASQUILL_D)
        @test_throws DomainError wind_speed(4.0, -1, SURFACE_WATER, PASQUILL_D)
        @test_throws DomainError wind_speed(-1.0, 50, SURFACE_WATER, PASQUILL_D)
    end

    @testset "lateral dispersion" begin
        for k in PASQUILL_CLASSES
            @test lateral_dispersion(0, k) == 0
            @test lateral_dispersion(100, k) > 0
            @test lateral_dispersion(100, k) < lateral_dispersion(1000, k)
            # No broadening at or below the reference duration.
            @test lateral_dispersion(500, k; release_duration = 60) ==
                  lateral_dispersion(500, k)
            @test lateral_dispersion(500, k; release_duration = 600) ==
                  lateral_dispersion(500, k)
            # Broadening above it, by exactly the stated factor.
            @test lateral_dispersion(500, k; release_duration = 3600) ≈
                  lateral_dispersion(500, k) * 6.0^0.2
        end
        # Unstable air disperses laterally faster than stable air.
        σ = [lateral_dispersion(1000, k) for k in PASQUILL_CLASSES]
        @test issorted(σ; rev = true)
        @test_throws DomainError lateral_dispersion(-1, PASQUILL_D)
        @test_throws DomainError lateral_dispersion(100, PASQUILL_D; release_duration = 0)
    end

    @testset "vertical dispersion" begin
        for k in PASQUILL_CLASSES, r in ROUGHNESS_CLASSES

            for x in (10.0, 100.0, 1000.0, 10_000.0)
                @test vertical_dispersion(x, k, r) > 0
            end
            # Monotone in distance over the range of interest.
            xs = [10.0, 100.0, 1000.0, 10_000.0]
            @test issorted([vertical_dispersion(x, k, r) for x in xs])
        end
        # Unstable air disperses vertically faster than stable air.
        σ = [vertical_dispersion(1000, k, ROUGHNESS_PASTURE) for k in PASQUILL_CLASSES]
        @test issorted(σ; rev = true)
        # Rougher ground mixes more.
        σ = [vertical_dispersion(1000, PASQUILL_D, r) for r in ROUGHNESS_CLASSES]
        @test issorted(σ)
        @test_throws DomainError vertical_dispersion(0, PASQUILL_D, ROUGHNESS_PASTURE)
        @test_throws DomainError vertical_dispersion(-1, PASQUILL_D, ROUGHNESS_PASTURE)
    end

    # Pasture is the neutral reference of the roughness correction: its
    # tabulated d₁, c₂ and d₂ all vanish, leaving F = ln(2.72) ≈ 1.
    @testset "pasture is the unit roughness correction" begin
        for x in (1.0, 100.0, 10_000.0)
            @test roughness_correction(x, ROUGHNESS_PASTURE) ≈ 1 atol = 1e-3
            @test vertical_dispersion(x, PASQUILL_D, ROUGHNESS_PASTURE) ≈
                  vertical_shape(x, PASQUILL_D) atol = 1e-3 * vertical_shape(x, PASQUILL_D)
        end
    end

    # The tables and formulae were transcribed from the 2021 implementation.
    # Reproduce it literally here and require exact agreement, so that a
    # mistyped coefficient cannot pass as a modelling choice.
    @testset "fidelity to the 2021 implementation" begin
        T1 = Dict(  # Tabel_1.csv
            'A' => (0.112, 1.06, 5.38e-4, 0.815),
            'B' => (0.13, 0.95, 6.52e-4, 0.75),
            'C' => (0.112, 0.92, 9.05e-4, 0.718),
            'D' => (0.098, 0.889, 1.35e-3, 0.688),
            'E' => (0.0609, 0.895, 1.96e-3, 0.684),
            'F' => (0.0638, 0.783, 1.36e-3, 0.672),
        )
        # Tabel_2.csv: z_0, c_1, d_1, c_2, d_2. The first two c_1 read 1.58
        # and 2.08 here, which is what the 2021 tables carried; Hosker
        # publishes 1.56 and 2.02. This reference implementation keeps the
        # 2021 values, so those two rows no longer match the package — see
        # the assertion below, which pins the size of the correction.
        T2 = [
            (0.01, 1.58, 0.048, 6.25e-4, 0.45),
            (0.04, 2.08, 0.0269, 7.76e-4, 0.37),
            (0.1, 2.72, 0.0, 0.0, 0.0),
            (0.4, 5.16, -0.098, 18.6, -0.225),
            (1.0, 7.37, -0.0957, 4.29e3, -0.6),
            (4.0, 11.7, -0.128, 4.59e4, -0.78),
        ]
        T3 = Dict(
            'A' => 0.22,
            'B' => 0.16,
            'C' => 0.11,
            'D' => 0.08,
            'E' => 0.06,
            'F' => 0.04,
        )
        T4 = Dict(  # Tabel_4.csv, by surface then class
            1 => Dict(
                'A' => 0.03,
                'B' => 0.05,
                'C' => 0.06,
                'D' => 0.08,
                'E' => 0.1,
                'F' => 0.12,
            ),
            2 => Dict(
                'A' => 0.1,
                'B' => 0.15,
                'C' => 0.2,
                'D' => 0.25,
                'E' => 0.35,
                'F' => 0.4,
            ),
            3 => Dict(
                'A' => 0.16,
                'B' => 0.24,
                'C' => 0.32,
                'D' => 0.4,
                'E' => 0.56,
                'F' => 0.64,
            ),
        )

        original_σ_z = function (x, P, s)
            a_1, b_1, a_2, b_2 = T1[P]
            g = (a_1 * x^b_1) / (1 + a_2 * x^b_2)
            z_0, c_1, d_1, c_2, d_2 = T2[s]
            F = if z_0 > 0.1
                log(c_1 * x^d_1 * (1 + (c_2 * x^d_2)^(-1)))
            else
                log((c_1 * x^d_1) / (1 + c_2 * x^d_2))
            end
            return g * F
        end
        original_σ_y = function (x, P, t_R)
            σy = (T3[P] * x) / (1 + 0.0001 * x)^0.5
            return t_R <= 600 ? σy : σy * (t_R / 600)^0.2
        end
        original_u_z = function (u_10, z, P, s)
            m = T4[s][P]
            return z <= 200 ? u_10 * (z / 10.0)^m : u_10 * (200 / 10.0)^m
        end

        for (i, k) in enumerate(PASQUILL_CLASSES)
            P = letter(k)
            for x in (1.0, 10.0, 137.0, 1000.0, 12_345.0, 100_000.0)
                for (j, r) in enumerate(ROUGHNESS_CLASSES)
                    if j ≤ 2
                        # Deliberate divergence, the only one in σ_z: the
                        # corrected coefficients give a smaller plume, by
                        # under 3 % over grassland and water and under 5 %
                        # over arable land, at every distance tested.
                        corrected = vertical_dispersion(x, k, r)
                        @test corrected < original_σ_z(x, P, j)
                        @test corrected / original_σ_z(x, P, j) > (j == 1 ? 0.97 : 0.95)
                    else
                        @test vertical_dispersion(x, k, r) == original_σ_z(x, P, j)
                    end
                end
                for t_R in (60.0, 600.0, 3600.0, 86_400.0)
                    @test lateral_dispersion(x, k; release_duration = t_R) ==
                          original_σ_y(x, P, t_R)
                end
            end
            for (j, s) in enumerate(WIND_PROFILE_SURFACES)
                for z in (1.0, 10.0, 50.3, 199.0, 200.0, 500.0)
                    @test wind_speed(4.0, z, s, k) == original_u_z(4.0, z, P, j)
                end
            end
        end
    end
end

@testset "Dispersion schemes" begin
    @testset "Briggs open country is the normative's σ_y" begin
        for k in PASQUILL_CLASSES, x in (0.0, 100.0, 1000.0, 10_000.0)

            @test briggs_lateral_dispersion(x, k, DISPERSION_BRIGGS_OPEN_COUNTRY) ==
                  lateral_dispersion(x, k)
            @test lateral_dispersion(x, k, DISPERSION_HOSKER) == lateral_dispersion(x, k)
        end
    end

    @testset "the scheme dispatch" begin
        for k in PASQUILL_CLASSES, r in ROUGHNESS_CLASSES, x in (100.0, 1000.0, 10_000.0)
            @test vertical_dispersion(x, k, DISPERSION_HOSKER, r) ==
                  vertical_dispersion(x, k, r)
            for s in (DISPERSION_BRIGGS_OPEN_COUNTRY, DISPERSION_BRIGGS_URBAN)
                @test vertical_dispersion(x, k, s, r) == briggs_vertical_dispersion(x, k, s)
                @test lateral_dispersion(x, k, s) == briggs_lateral_dispersion(x, k, s)
            end
            @test vertical_dispersion(x, k, DISPERSION_EIMUTIS_KONICEK, r) ==
                  eimutis_konicek_vertical_dispersion(x, k)
            @test lateral_dispersion(x, k, DISPERSION_EIMUTIS_KONICEK) ==
                  eimutis_konicek_lateral_dispersion(x, k)
            p = dispersion_parameters(
                x, k, DISPERSION_BRIGGS_URBAN, r; release_duration = 3600,)
            @test p.σy ==
                  meander_broadened(briggs_lateral_dispersion(x, k, DISPERSION_BRIGGS_URBAN), 3600)
            @test p.σz == briggs_vertical_dispersion(x, k, DISPERSION_BRIGGS_URBAN)
        end
        @test_throws ArgumentError briggs_lateral_coefficients(PASQUILL_D, DISPERSION_HOSKER)
        @test_throws ArgumentError briggs_vertical_dispersion(
            100.0, PASQUILL_D, DISPERSION_EIMUTIS_KONICEK,)
        @test_throws DomainError briggs_lateral_dispersion(-1.0, PASQUILL_D, DISPERSION_BRIGGS_URBAN)
        @test_throws DomainError briggs_vertical_dispersion(0.0, PASQUILL_D, DISPERSION_BRIGGS_URBAN)
    end

    @testset "meander broadening" begin
        @test meander_broadened(100.0, 60) == 100.0
        @test meander_broadened(100.0, 600) == 100.0
        @test meander_broadened(100.0, 3600) ≈ 100 * 6^0.2
        @test_throws DomainError meander_broadened(100.0, 0)
    end

    @testset "every scheme orders the classes and grows with distance" begin
        xs = (100.0, 300.0, 1000.0, 3000.0, 10_000.0)
        for s in DISPERSION_SCHEMES
            for k in PASQUILL_CLASSES
                @test issorted([lateral_dispersion(x, k, s) for x in xs])
                @test issorted([vertical_dispersion(x, k, s, ROUGHNESS_PASTURE) for x in xs])
            end
            # Weakly ordered: the urban set pairs A with B and E with F.
            for x in xs
                @test issorted([lateral_dispersion(x, k, s) for k in PASQUILL_CLASSES]; rev = true)
                @test issorted(
                    [vertical_dispersion(x, k, s, ROUGHNESS_PASTURE)
                     for k in PASQUILL_CLASSES];
                    rev = true,
                )
            end
        end
    end

    # The urban σ_z exceeds the open-country one at every distance of the
    # band; the urban σ_y does so only in the near field, its `(1 + 4×10⁻⁴x)`
    # denominator growing four times faster than the open-country one.
    @testset "urban air disperses faster than open country" begin
        for k in PASQUILL_CLASSES
            for x in (100.0, 1000.0, 10_000.0)
                @test briggs_vertical_dispersion(x, k, DISPERSION_BRIGGS_URBAN) >
                      briggs_vertical_dispersion(x, k, DISPERSION_BRIGGS_OPEN_COUNTRY)
            end
            for x in (100.0, 1000.0)
                @test briggs_lateral_dispersion(x, k, DISPERSION_BRIGGS_URBAN) >
                      briggs_lateral_dispersion(x, k, DISPERSION_BRIGGS_OPEN_COUNTRY)
            end
        end
    end

    @testset "the Eimutis–Konicek pieces join" begin
        for k in PASQUILL_CLASSES, x in EIMUTIS_KONICEK_RANGES

            below = eimutis_konicek_vertical_dispersion(prevfloat(x), k)
            above = eimutis_konicek_vertical_dispersion(x, k)
            @test isapprox(below, above; rtol = 0.02)
        end
        @test eimutis_konicek_vertical_coefficients(PASQUILL_A, 0.0) ===
              eimutis_konicek_vertical_coefficients(PASQUILL_A, 99.0)
        @test_throws DomainError eimutis_konicek_vertical_dispersion(0.0, PASQUILL_D)
        @test_throws DomainError eimutis_konicek_vertical_coefficients(PASQUILL_D, -1.0)
        @test eimutis_konicek_lateral_dispersion(0.0, PASQUILL_D) == 0
        @test_throws DomainError eimutis_konicek_lateral_dispersion(-1.0, PASQUILL_D)
    end

    @testset "validity ranges" begin
        for s in DISPERSION_SCHEMES
            lo, hi = validity_range(s)
            @test lo == 100
            @test within_validity(lo, s) && within_validity(hi, s)
            @test !within_validity(prevfloat(lo), s) && !within_validity(nextfloat(hi), s)
        end
        @test validity_range(DISPERSION_HOSKER)[2] == 10_000
        @test validity_range(DISPERSION_BRIGGS_OPEN_COUNTRY)[2] == 10_000
        @test validity_range(DISPERSION_BRIGGS_URBAN)[2] == 10_000
        @test validity_range(DISPERSION_EIMUTIS_KONICEK)[2] == 100_000
    end

    @testset "the layer mean of the wind profile" begin
        u₁₀ = 4.0
        for s in WIND_PROFILE_SURFACES, k in PASQUILL_CLASSES

            for (z₁, z₂) in ((10.0, 60.0), (50.0, 250.0), (210.0, 400.0), (5.0, 15.0))
                q, _ = quadgk(z -> wind_speed(u₁₀, z, s, k), z₁, z₂)
                @test layer_mean_wind_speed(u₁₀, z₁, z₂, s, k) ≈ q / (z₂ - z₁) rtol = 1e-10
            end
            @test layer_mean_wind_speed(u₁₀, 50.0, 50.0, s, k) ==
                  wind_speed(u₁₀, 50.0, s, k)
            # Positive shear: the mean over a layer lies between the winds at its ends.
            @test wind_speed(u₁₀, 50.0, s, k) ≤
                  layer_mean_wind_speed(u₁₀, 50.0, 100.0, s, k) ≤
                  wind_speed(u₁₀, 100.0, s, k)
        end
        @test_throws DomainError layer_mean_wind_speed(
            u₁₀, 0.0, 10.0, SURFACE_WATER, PASQUILL_D,)
        @test_throws DomainError layer_mean_wind_speed(
            u₁₀, 20.0, 10.0, SURFACE_WATER, PASQUILL_D,)
        @test_throws DomainError layer_mean_wind_speed(
            -1.0, 10.0, 20.0, SURFACE_WATER, PASQUILL_D,)
    end
end
