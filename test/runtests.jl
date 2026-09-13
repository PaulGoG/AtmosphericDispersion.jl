using AtmosphericDispersion
using Test
using Aqua
using ExplicitImports
using JET
using QuadGK
using TOML

@testset "AtmosphericDispersion.jl" begin
    @testset "Static quality assurance" begin
        Aqua.test_all(AtmosphericDispersion)
        @test check_no_implicit_imports(AtmosphericDispersion) === nothing
        @test check_no_stale_explicit_imports(AtmosphericDispersion) === nothing
        @test check_all_qualified_accesses_via_owners(AtmosphericDispersion) === nothing
        # Scoped to this package: unscoped, JET walks into QuadGK's generic
        # Gauss-Kronrod machinery and reports unreachable branches of
        # LinearAlgebra's norm, which say nothing about the code here.
        JET.test_package(AtmosphericDispersion; target_modules = (AtmosphericDispersion,))
    end

    # Guard the boundary that report is really about: the quadrature call must
    # stay concretely typed, or the whole of QuadGK resolves dynamically.
    @testset "Inference at the quadrature boundary" begin
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
        @test Base.return_types(depletion_integral, (Float64, Site, PasquillClass)) ==
              [Float64]
        @test @inferred(depletion_integral(1000.0, site, PASQUILL_D)) isa Float64
        @test @inferred(dry_depletion_factor(1000.0, site, PASQUILL_D, TRITIATED_WATER)) isa
              Float64
        @test @inferred(effective_height(1000.0, site, PASQUILL_D)) isa Float64
        @test @inferred(corrected_vertical_dispersion(1000.0, site, PASQUILL_D)) isa Float64
    end

    @testset "SectorGrid" begin
        @testset "construction" begin
            @test nsectors(SectorGrid(16)) == 16
            @test sector_width(SectorGrid(16)) ≈ deg2rad(22.5)
            @test sector_width(SectorGrid(8)) ≈ deg2rad(45.0)
            @test_throws ArgumentError SectorGrid(2)
            @test_throws ArgumentError SectorGrid(15)
        end

        @testset "bearing" begin
            @test bearing(0, 1) ≈ 0                # north
            @test bearing(1, 0) ≈ π / 2            # east
            @test bearing(0, -1) ≈ π               # south
            @test bearing(-1, 0) ≈ 3π / 2          # west
            @test bearing(1, 1) ≈ π / 4            # north-east
            @test bearing(0, 0) == 0
            @test 0 ≤ bearing(-1, -1) < 2π
        end

        # The 2021 implementation binned by truncating a ratio, which put every
        # cardinal direction on a sector edge and let one ulp of round-off decide
        # the bin: ESE and SE collapsed into one sector, SSW and SW into another,
        # and two sectors became unreachable. Each cardinal direction must land in
        # its own sector, at that sector's centre.
        @testset "the sixteen cardinal directions" begin
            g = SectorGrid(16)
            indices = [sector_of(g, deg2rad(22.5 * (i - 1))) for i = 1:16]

            @test indices == collect(1:16)
            @test length(unique(indices)) == 16
            @test [sector_name(g, k) for k in indices] == collect(CARDINAL_16)

            for k = 1:16
                @test sector_bearing(g, k) ≈ deg2rad(22.5 * (k - 1))
                @test sector_of(g, sector_bearing(g, k)) == k
            end
        end

        @testset "cardinal directions from Cartesian displacements" begin
            g = SectorGrid(16)
            # x east, y north; a receptor due east of the source is in sector 5.
            @test sector_of(g, 0.0, 1.0) == 1     # N
            @test sector_of(g, 1.0, 1.0) == 3     # NE
            @test sector_of(g, 1.0, 0.0) == 5     # E
            @test sector_of(g, 1.0, -1.0) == 7    # SE
            @test sector_of(g, 0.0, -1.0) == 9    # S
            @test sector_of(g, -1.0, -1.0) == 11  # SW
            @test sector_of(g, -1.0, 0.0) == 13   # W
            @test sector_of(g, -1.0, 1.0) == 15   # NW
            # Distance must not affect the sector.
            @test sector_of(g, 1e-9, 0.0) == sector_of(g, 1e9, 0.0)
        end

        @testset "sector boundaries" begin
            g = SectorGrid(16)
            half = sector_width(g) / 2
            for k = 1:16
                c = sector_bearing(g, k)
                lo, hi = sector_bounds(g, k)
                @test mod2pi(c - half) ≈ lo
                @test mod2pi(c + half) ≈ hi
                # Interior points either side of centre stay in the sector.
                @test sector_of(g, mod2pi(c - 0.99half)) == k
                @test sector_of(g, mod2pi(c + 0.99half)) == k
                # A bearing exactly on a boundary belongs to one of the two
                # sectors meeting there. The tie rule takes the higher, but a
                # boundary bearing is generally not representable, so round-off
                # may place it either side. Both answers are defensible.
                @test sector_of(g, mod2pi(c + half)) in (k, mod(k, 16) + 1)
            end
        end

        @testset "wrapping" begin
            g = SectorGrid(16)
            # Away from the boundaries the sector is invariant under whole turns.
            for k = 1:16, offset in (-0.4, -0.1, 0.0, 0.1, 0.4)
                β = sector_bearing(g, k) + offset * sector_width(g)
                @test sector_of(g, β) == k
                @test sector_of(g, β + 2π) == k
                @test sector_of(g, β - 2π) == k
                @test sector_of(g, β + 6π) == k
            end
            # Any bearing whatsoever lands in a valid sector.
            for β in range(-4π, 4π; length = 401)
                @test 1 ≤ sector_of(g, β) ≤ 16
            end
        end

        @testset "the sectors partition the circle" begin
            g = SectorGrid(16)
            # Sample midway between boundaries so no sample sits on one.
            counts = zeros(Int, 16)
            m = 16 * 64
            for j = 0:(m-1)
                counts[sector_of(g, 2π * (j + 0.5) / m)] += 1
            end
            @test all(==(64), counts)
        end

        @testset "opposite" begin
            for n in (4, 8, 16, 36)
                g = SectorGrid(n)
                for k = 1:n
                    o = opposite(g, k)
                    @test 1 ≤ o ≤ n
                    @test opposite(g, o) == k
                    @test o != k
                    # Opposite sectors are half a turn apart.
                    δ = abs(sector_bearing(g, k) - sector_bearing(g, o))
                    @test min(δ, 2π - δ) ≈ π
                end
            end
            g = SectorGrid(16)
            @test sector_name(g, opposite(g, 1)) == "S"
            @test sector_name(g, opposite(g, 5)) == "W"
            @test_throws BoundsError opposite(g, 0)
            @test_throws BoundsError opposite(g, 17)
        end

        @testset "naming" begin
            g = SectorGrid(16)
            @test sector_name(g, 1) == "N"
            @test sector_name(g, 5) == "E"
            @test sector_name(g, 9) == "S"
            @test sector_name(g, 13) == "W"
            @test sector_name(SectorGrid(8), 3) == "90.0°"
        end
    end

    @testset "PasquillClass" begin
        @test length(PASQUILL_CLASSES) == 6
        @test PASQUILL_CLASSES[1] === PASQUILL_A
        @test PASQUILL_CLASSES[6] === PASQUILL_F
        @test pasquill('D') === PASQUILL_D
        @test pasquill('d') === PASQUILL_D
        @test pasquill("F") === PASQUILL_F
        @test pasquill(PASQUILL_B) === PASQUILL_B
        @test_throws ArgumentError pasquill('G')
        @test_throws ArgumentError pasquill("DE")
        for (i, c) in enumerate(PASQUILL_CLASSES)
            @test classindex(c) == i
            @test pasquill(letter(c)) === c
        end
    end

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
            @test_throws DomainError lateral_dispersion(
                100,
                PASQUILL_D;
                release_duration = 0,
            )
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
                      vertical_shape(x, PASQUILL_D) atol =
                    1e-3 * vertical_shape(x, PASQUILL_D)
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
            @test_throws DomainError final_buoyant_rise(
                buoyancy_flux(dense, air),
                5.0,
                1e-3,
            )
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

    @testset "Plume rise" begin
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
            neutral = Atmosphere(;
                reference_speed = 4.0,
                temperature = 287.0,
                density = 1.2,
                lapse_rate = -STANDARD_GRAVITY / DRY_AIR_SPECIFIC_HEAT,
                surface = SURFACE_AGRICULTURAL,
                roughness = ROUGHNESS_PASTURE,
            )
            @test final_buoyant_rise(F, 4.0, stability_parameter(neutral)) ≥
                  final_buoyant_rise(F, 4.0, S)
            # Unstable air has no stable ceiling, and must not raise a domain
            # error from a fractional power of a negative stability parameter.
            unstable = Atmosphere(;
                reference_speed = 4.0,
                temperature = 287.0,
                density = 1.2,
                lapse_rate = -0.02,
                surface = SURFACE_AGRICULTURAL,
                roughness = ROUGHNESS_PASTURE,
            )
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
                (x < 3.5 * original_X_0(F) && hbtranzitie <= hbfinal) ? hbtranzitie :
                hbfinal
            end
            original_hm_final = function (Fm, w_0, D, u, S)
                min(
                    1.5 * w_0 * D / u,
                    4 * (Fm / S)^(1 / 4),
                    1.5 * (Fm / u)^(1 / 3) * S^(-1 / 6),
                )
            end
            original_hm = function (x, Fm, w_0, D, u, S)
                hmfinal = original_hm_final(Fm, w_0, D, u, S)
                hmtranzitie = 1.89 * (w_0^2 * D / (u * (w_0 + 3u)))^(2 / 3) * x^(1 / 3)
                hmtranzitie <= hmfinal ? hmtranzitie : hmfinal
            end
            original_hmb = function (x, F, Fm, w_0, D, u, S)
                hmbfinal = original_hm_final(Fm, w_0, D, u, S) + original_hb_final(F, u, S)
                hmbtranzitie =
                    3^(1 / 3) *
                    (Fm * x / ((1 / 3 + u / w_0)^2 * u^2) + F * x^2 / (0.5 * u^3))^(1 / 3)
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
                @test final_buoyant_rise(F, u, Sv, THESIS_RISE) ==
                      original_hb_final(F, u, Sv)
                @test final_momentum_rise(Fₘ, w₀, D, u, Sv, THESIS_RISE) ==
                      original_hm_final(Fₘ, w₀, D, u, Sv)
                for x in (1.0, 10.0, 137.0, 1000.0, 10_000.0)
                    @test buoyant_rise(x, F, u, Sv, THESIS_RISE) == original_hb(x, F, u, Sv)
                    @test momentum_rise(x, Fₘ, w₀, D, u, Sv, THESIS_RISE) ==
                          original_hm(x, Fₘ, w₀, D, u, Sv)
                    @test combined_rise(x, F, Fₘ, w₀, u, Sv, D, THESIS_RISE) ==
                          original_hmb(x, F, Fₘ, w₀, D, u, Sv)
                end
            end
        end
    end

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

    @testset "Dilution" begin
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

        @testset "plume frame" begin
            # Wind from the north transports to the south.
            x, y = plume_frame(0.0, -1000.0, 0.0)
            @test x ≈ 1000 && abs(y) < 1e-9
            @test plume_frame(0.0, 1000.0, 0.0)[1] ≈ -1000        # upwind
            # Wind from the west transports to the east.
            x, y = plume_frame(1000.0, 0.0, 3π / 2)
            @test x ≈ 1000 && abs(y) < 1e-9
            # The frame is orthonormal: it preserves distance.
            for β in range(0, 2π; length = 17), (e, n) in ((300.0, 400.0), (-120.0, 50.0))
                x, y = plume_frame(e, n, β)
                @test hypot(x, y) ≈ hypot(e, n)
            end
        end

        @testset "instantaneous field" begin
            # Nothing upwind.
            @test dilution_instantaneous(0.0, 1000.0, 0.0, site, PASQUILL_D, 0.0) == 0
            # Peaked on the axis and symmetric across it.
            onaxis = dilution_instantaneous(0.0, -1000.0, 0.0, site, PASQUILL_D, 0.0)
            @test onaxis > 0
            for off in (50.0, 200.0)
                left = dilution_instantaneous(off, -1000.0, 0.0, site, PASQUILL_D, 0.0)
                right = dilution_instantaneous(-off, -1000.0, 0.0, site, PASQUILL_D, 0.0)
                @test left ≈ right
                @test left < onaxis
            end
            # Falls off downwind, far from the source.
            far = [
                dilution_instantaneous(0.0, -r, 0.0, site, PASQUILL_D, 0.0) for
                r in (2_000.0, 5_000.0, 20_000.0)
            ]
            @test issorted(far; rev = true)
            # Equivariant under a common rotation of receptor and wind.
            for β in range(0, 2π; length = 13)
                rotated = dilution_instantaneous(
                    1000sin(β + π),
                    1000cos(β + π),
                    0.0,
                    site,
                    PASQUILL_D,
                    β,
                )
                @test rotated ≈ onaxis
            end
            @test_throws DomainError dilution_instantaneous(
                0.0,
                -1000.0,
                -1.0,
                site,
                PASQUILL_D,
                0.0,
            )
        end

        @testset "extended field" begin
            # Uniform across the sector, zero outside it.
            onaxis = dilution_extended(0.0, -1000.0, site, PASQUILL_D, 0.0)
            @test onaxis > 0
            inside = dilution_extended(90.0, -1000.0, site, PASQUILL_D, 0.0)
            @test inside ≈ onaxis        # crosswind-uniform within the sector
            @test dilution_extended(1000.0, 0.0, site, PASQUILL_D, 0.0) == 0   # 90° off
            @test dilution_extended(0.0, 1000.0, site, PASQUILL_D, 0.0) == 0   # upwind
            # A wider sector spreads the same material further, so dilutes more.
            @test dilution_extended(0.0, -1000.0, site, PASQUILL_D, 0.0, SectorGrid(8)) <
                  onaxis
        end

        @testset "long-term field" begin
            g = SectorGrid(16)
            stab = [0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158]
            stab ./= sum(stab)

            # A uniform rose gives a rotationally symmetric field.
            uniform = WindRose(g, fill(1 / 16, 16), BlowingToward(); stability = stab)
            values = [
                dilution_long_term(
                    1000sin(sector_bearing(g, k)),
                    1000cos(sector_bearing(g, k)),
                    site,
                    uniform,
                ) for k = 1:16
            ]
            @test all(v -> v ≈ first(values), values)
            @test first(values) > 0

            # The two conventions give fields that are exact mirror images,
            # sector by sector. This is the defect, isolated.
            f = [
                0.07849,
                0.07849,
                0.07221,
                0.07378,
                0.07849,
                0.05965,
                0.05651,
                0.0471,
                0.05181,
                0.05024,
                0.04867,
                0.04553,
                0.05495,
                0.06436,
                0.06279,
                0.07692,
            ]
            f ./= sum(f)
            from = WindRose(g, f, BlowingFrom(); stability = stab)
            toward = WindRose(g, f, BlowingToward(); stability = stab)
            χ(rose, k) = dilution_long_term(
                1000sin(sector_bearing(g, k)),
                1000cos(sector_bearing(g, k)),
                site,
                rose,
            )
            for k = 1:16
                @test χ(from, k) ≈ χ(toward, opposite(g, k))
            end
            # And they disagree substantially where the rose is asymmetric.
            ratios = [χ(toward, k) / χ(from, k) for k = 1:16]
            @test maximum(ratios) > 1.4
            @test minimum(ratios) < 0.7
            # The most exposed sector is exactly opposite between the two.
            @test argmax([χ(toward, k) for k = 1:16]) == opposite(g, argmax([χ(from, k) for k = 1:16]))

            @test dilution_long_term(0.0, 0.0, site, toward) == 0
            # A sector the wind never blows towards receives nothing.
            single = zeros(16)
            single[1] = 1.0
            only_north = WindRose(g, single, BlowingToward(); stability = stab)
            @test χ(only_north, 1) > 0
            @test χ(only_north, 9) == 0
        end
    end

    @testset "Depletion" begin
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
        inert = Nuclide(;
            name = "inert",
            decay_constant = 0.0,
            deposition_velocity = DepositionVelocity(0.0, 0.0),
        )

        @testset "nuclides" begin
            @test TRITIATED_WATER.decay_constant == TRITIUM_DECAY_CONSTANT
            @test TRITIUM_GAS.decay_constant == TRITIUM_DECAY_CONSTANT
            # HTO deposits an order of magnitude more readily than HT: exactly
            # a factor of ten in the low velocity, sixteen in the high.
            @test TRITIATED_WATER.deposition_velocity.low ≈
                  10 * TRITIUM_GAS.deposition_velocity.low
            @test TRITIATED_WATER.deposition_velocity.high ≈
                  16 * TRITIUM_GAS.deposition_velocity.high
            # Tritium's half-life is about 12.3 years.
            @test log(2) / TRITIUM_DECAY_CONSTANT / (365.25 * 86400) ≈ 12.3 atol = 0.05
            @test_throws ArgumentError DepositionVelocity(0.5, 0.1)
            @test_throws ArgumentError DepositionVelocity(-0.1, 0.1)
            @test_throws ArgumentError Nuclide(;
                name = "x",
                decay_constant = -1.0,
                deposition_velocity = DepositionVelocity(0.0, 0.0),
            )
            @test occursin("HTO", sprint(show, TRITIATED_WATER))
        end

        @testset "decay" begin
            u = transport_wind_speed(site, PASQUILL_D)
            @test decay_factor(0.0, u, TRITIATED_WATER) == 1
            @test decay_factor(1000.0, u, TRITIATED_WATER) ==
                  exp(-TRITIUM_DECAY_CONSTANT * 1000 / u)
            @test decay_factor(1000.0, u, inert) == 1
            # Tritium is far too long-lived for transit decay to matter.
            @test decay_factor(1e5, u, TRITIATED_WATER) > 0.9999
            @test issorted(
                [decay_factor(x, u, TRITIATED_WATER) for x in (0.0, 1e4, 1e6, 1e8)];
                rev = true,
            )
            @test_throws DomainError decay_factor(-1.0, u, TRITIATED_WATER)
            @test_throws DomainError decay_factor(100.0, 0.0, TRITIATED_WATER)
        end

        @testset "dry depletion" begin
            for class in PASQUILL_CLASSES
                @test depletion_integral(DEPLETION_INTEGRAL_FLOOR, site, class) == 0
                @test depletion_integral(0.5, site, class) == 0
                # The integral accumulates with distance.
                I = [depletion_integral(x, site, class) for x in (10.0, 1e3, 1e4, 1e5)]
                @test issorted(I)
                @test all(≥(0), I)
                # A surviving fraction, so in (0, 1] and falling with distance.
                f = [
                    dry_depletion_factor(x, site, class, TRITIATED_WATER) for
                    x in (10.0, 1e3, 1e4, 1e5)
                ]
                @test all(v -> 0 < v ≤ 1, f)
                @test issorted(f; rev = true)
                # No deposition velocity, no depletion.
                @test dry_depletion_factor(1e4, site, class, inert) == 1
                # HT deposits less than HTO, so more of it survives.
                @test dry_depletion_factor(1e4, site, class, TRITIUM_GAS) >
                      dry_depletion_factor(1e4, site, class, TRITIATED_WATER)
            end
        end

        @testset "wet depletion" begin
            @test wet_depletion_factor(0.0, PRECIPITATION_RAIN, 1.0) == 1
            for p in PRECIPITATION_TYPES, r in PRECIPITATION_RATES
                f = [wet_depletion_factor(t, p, r) for t in (0.0, 600.0, 3600.0, 86400.0)]
                @test all(v -> 0 < v ≤ 1, f)
                @test issorted(f; rev = true)
            end
            # Heavier rain washes out faster; snow barely at all.
            @test wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 5.0) <
                  wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 0.5)
            @test wet_depletion_factor(3600.0, PRECIPITATION_SNOW, 5.0) >
                  wet_depletion_factor(3600.0, PRECIPITATION_RAIN, 0.5)
            @test_throws DomainError wet_depletion_factor(-1.0, PRECIPITATION_RAIN, 1.0)
        end

        # The composition is multiplicative, not additive. With nothing removed
        # at all every factor is one and so is the product. The 2021 code added
        # the wet and dry factors and returned two in this limit.
        @testset "the no-depletion limit is one" begin
            for class in PASQUILL_CLASSES, x in (10.0, 1e3, 1e5)
                @test depletion_factor(x, site, class, inert) == 1
                @test depletion_factor(x, site, class, inert; washout_duration = 0.0) == 1
            end
        end

        @testset "composition" begin
            for class in PASQUILL_CLASSES
                x, t = 1e4, 3600.0
                u = transport_wind_speed(site, class)
                expected =
                    decay_factor(x, u, TRITIATED_WATER) *
                    dry_depletion_factor(x, site, class, TRITIATED_WATER) *
                    wet_depletion_factor(t, PRECIPITATION_RAIN, 1.0)
                @test depletion_factor(
                    x,
                    site,
                    class,
                    TRITIATED_WATER;
                    washout_duration = t,
                    precipitation = PRECIPITATION_RAIN,
                    rate = 1.0,
                ) ≈ expected
                @test 0 < expected ≤ 1
                # Rain can only remove more.
                @test depletion_factor(
                    x,
                    site,
                    class,
                    TRITIATED_WATER;
                    washout_duration = t,
                ) < depletion_factor(x, site, class, TRITIATED_WATER)
            end
        end

        @testset "depletion inside the long-term class sum" begin
            g = SectorGrid(16)
            stab = [0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158]
            stab ./= sum(stab)
            rose = WindRose(g, fill(1 / 16, 16), BlowingToward(); stability = stab)
            bare = dilution_long_term(0.0, -5000.0, site, rose)
            # An inert species leaves the field untouched: the depletion factor
            # is exactly one, not six as the unweighted sum of the 2021 code
            # would have given.
            @test dilution_long_term(0.0, -5000.0, site, rose; nuclide = inert) == bare
            depleted =
                dilution_long_term(0.0, -5000.0, site, rose; nuclide = TRITIATED_WATER)
            @test 0 < depleted < bare
            @test depleted / bare > 0.9        # tritium depletes slowly
            # Washout removes more still.
            @test dilution_long_term(
                0.0,
                -5000.0,
                site,
                rose;
                nuclide = TRITIATED_WATER,
                washout_duration = 3600.0,
            ) < depleted
        end
    end

    @testset "Deposition and resuspension" begin
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
        Q = 1e12

        @testset "dry deposition" begin
            @test dry_deposition(0.0, TRITIATED_WATER) == 0
            @test dry_deposition(1000.0, TRITIATED_WATER) ==
                  TRITIATED_WATER.deposition_velocity.high * 1000
            # Deposition takes the high velocity, depletion the low one.
            @test dry_deposition(1000.0, TRITIATED_WATER) >
                  TRITIATED_WATER.deposition_velocity.low * 1000
            # HT deposits far less than HTO.
            @test dry_deposition(1000.0, TRITIUM_GAS) <
                  dry_deposition(1000.0, TRITIATED_WATER)
            @test_throws DomainError dry_deposition(-1.0, TRITIATED_WATER)
        end

        @testset "wet deposition" begin
            ω = wet_deposition(
                0.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
            )
            @test ω > 0
            # Nothing upwind, nothing for a null release.
            @test wet_deposition(
                0.0,
                1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
            ) == 0
            @test wet_deposition(
                0.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = 0.0,
            ) == 0
            # Linear in released activity.
            @test wet_deposition(
                0.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = 2Q,
                rate = 1.0,
            ) ≈ 2ω
            # Peaked on the axis, symmetric across it.
            off = wet_deposition(
                300.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
            )
            @test off < ω
            @test off ≈ wet_deposition(
                -300.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
            )
            # Heavier rain deposits more; snow far less.
            @test wet_deposition(
                0.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
                rate = 5.0,
            ) > ω
            @test wet_deposition(
                0.0,
                -1000.0,
                site,
                PASQUILL_D,
                0.0,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
                precipitation = PRECIPITATION_SNOW,
            ) < ω / 100

            # The column normalisation is √(2π) Σ_y u, which is what integrating
            # the Gaussian plume over all z leaves. The 2021 code wrote √2 π,
            # larger by exactly √π, understating deposition by that factor.
            Λ = washout_coefficients(PRECIPITATION_RAIN, 1.0).high
            u = transport_wind_speed(site, PASQUILL_D)
            x, _ = plume_frame(0.0, -1000.0, 0.0)
            Σy = corrected_lateral_dispersion(x, site, PASQUILL_D)
            expected = Λ * Q * decay_factor(x, u, TRITIATED_WATER) / (sqrt(2π) * Σy * u)
            @test ω ≈ expected
            @test sqrt(2) * π / sqrt(2π) ≈ sqrt(π)
            @test ω ≈
                  sqrt(π) * Λ * Q * decay_factor(x, u, TRITIATED_WATER) /
                  (sqrt(2) * π * Σy * u)
        end

        @testset "sector-averaged wet deposition" begin
            ω = wet_deposition_sector(
                1000.0,
                site,
                PASQUILL_D,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
            )
            @test ω > 0
            @test wet_deposition_sector(
                0.0,
                site,
                PASQUILL_D,
                TRITIATED_WATER;
                activity = Q,
            ) == 0
            # Spread over an arc that grows with distance, so falls as 1/r.
            @test wet_deposition_sector(
                2000.0,
                site,
                PASQUILL_D,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
            ) ≈ ω / 2 rtol = 1e-6
            # A wider sector spreads the same material further.
            @test wet_deposition_sector(
                1000.0,
                site,
                PASQUILL_D,
                TRITIATED_WATER;
                activity = Q,
                rate = 1.0,
                sectors = SectorGrid(8),
            ) ≈ ω / 2 rtol = 1e-6
        end

        @testset "resuspension" begin
            c = RESUSPENSION_COEFFICIENTS
            @test resuspension_factor(0.0) ≈ c.A + c.B
            # Falls monotonically, and stays positive.
            ks = [resuspension_factor(t) for t in (0.0, 1.0, 30.0, 365.0, 3650.0)]
            @test issorted(ks; rev = true)
            @test all(>(0), ks)
            # About a factor 38 over the first year, four orders over ten.
            @test resuspension_factor(0.0) / resuspension_factor(365.0) ≈ 38 atol = 1
            @test resuspension_factor(0.0) / resuspension_factor(3650.0) > 1e4
            # The slow term is all that survives in the long run.
            @test resuspension_factor(1e5) ≈ c.B * exp(-c.λ₂ * 1e5)
            @test_throws DomainError resuspension_factor(-1.0)

            @test resuspended_concentration(0.0, 10.0) == 0
            @test resuspended_concentration(1000.0, 0.0) ≈ 1000 * (c.A + c.B)
            @test_throws DomainError resuspended_concentration(-1.0, 10.0)
        end

        @testset "depletion reaches the short-range regimes too" begin
            inert = Nuclide(;
                name = "inert",
                decay_constant = 0.0,
                deposition_velocity = DepositionVelocity(0.0, 0.0),
            )
            bare = dilution_instantaneous(0.0, -5000.0, 0.0, site, PASQUILL_D, 0.0)
            @test dilution_instantaneous(
                0.0,
                -5000.0,
                0.0,
                site,
                PASQUILL_D,
                0.0;
                nuclide = inert,
            ) == bare
            @test 0 <
                  dilution_instantaneous(
                      0.0,
                      -5000.0,
                      0.0,
                      site,
                      PASQUILL_D,
                      0.0;
                      nuclide = TRITIATED_WATER,
                  ) <
                  bare
            bare_e = dilution_extended(0.0, -5000.0, site, PASQUILL_D, 0.0)
            @test dilution_extended(0.0, -5000.0, site, PASQUILL_D, 0.0; nuclide = inert) ==
                  bare_e
            @test 0 <
                  dilution_extended(
                      0.0,
                      -5000.0,
                      site,
                      PASQUILL_D,
                      0.0;
                      nuclide = TRITIATED_WATER,
                  ) <
                  bare_e
        end
    end

    @testset "Configuration" begin
        reference = joinpath(@__DIR__, "..", "config", "reference.toml")

        @testset "the reference configuration loads" begin
            @test isfile(reference)
            c = load_configuration(reference)
            @test c isa RunConfiguration
            @test c.site.source.height == 50.3
            @test c.site.atmosphere.surface === SURFACE_AGRICULTURAL
            @test c.site.atmosphere.roughness === ROUGHNESS_PASTURE
            @test c.nuclide.name == "HTO"
            @test c.nuclide.decay_constant == TRITIUM_DECAY_CONSTANT
            @test nsectors(grid(c.rose)) == 16
            @test c.activity == 6.8e14
            @test c.precipitation === PRECIPITATION_RAIN
            @test c.spacing < c.extent
            # It drives the solver.
            @test dilution_long_term(0.0, -c.extent, c.site, c.rose; nuclide = c.nuclide) >
                  0
            @test_throws ArgumentError load_configuration(reference * ".missing")
        end

        # Mutate a copy of the reference table and require the named failure.
        base() = TOML.parsefile(reference)
        function withkey(f)
            t = base()
            f(t)
            return t
        end

        @testset "missing keys are named" begin
            for (table, key) in (
                ("source", "height"),
                ("source", "diameter"),
                ("atmosphere", "temperature"),
                ("atmosphere", "surface"),
                ("nuclide", "decay_constant"),
                ("release", "activity"),
                ("grid", "extent"),
                ("wind_rose", "convention"),
            )
                t = withkey(t -> delete!(t[table], key))
                err = try
                    configuration_from(t)
                    nothing
                catch e
                    e
                end
                @test err isa ConfigurationError
                @test err.path == "$table.$key"
                @test occursin("missing", err.message)
            end
            # A missing table is named too.
            t = withkey(t -> delete!(t, "release"))
            err = try
                configuration_from(t)
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "release"
        end

        @testset "enumerated choices list their options" begin
            for (table, key, options) in (
                ("atmosphere", "surface", "agricultural"),
                ("atmosphere", "roughness", "pasture"),
                ("wind_rose", "convention", "blowing_from"),
                ("precipitation", "type", "rain"),
            )
                t = withkey(t -> (t[table][key] = "nonsense"))
                err = try
                    configuration_from(t)
                    nothing
                catch e
                    e
                end
                @test err isa ConfigurationError
                @test err.path == "$table.$key"
                @test occursin("must be one of", err.message)
                @test occursin(options, err.message)
            end
        end

        @testset "numerical bounds" begin
            cases = [
                (t -> (t["source"]["height"] = 0.0), "source.height"),
                (t -> (t["source"]["diameter"] = -1.0), "source.diameter"),
                (t -> (t["source"]["exit_velocity"] = -1.0), "source.exit_velocity"),
                (t -> (t["atmosphere"]["temperature"] = 0.0), "atmosphere.temperature"),
                (t -> (t["atmosphere"]["density"] = -1.0), "atmosphere.density"),
                (t -> (t["release"]["activity"] = -1.0), "release.activity"),
                (t -> (t["release"]["duration"] = 0.0), "release.duration"),
                (t -> (t["grid"]["extent"] = -1.0), "grid.extent"),
                (t -> (t["grid"]["spacing"] = 0.0), "grid.spacing"),
                (
                    t -> (t["buildings"]["wake_coefficient"] = -1.0),
                    "buildings.wake_coefficient",
                ),
                (
                    t -> (t["precipitation"]["washout_duration"] = -1.0),
                    "precipitation.washout_duration",
                ),
            ]
            for (mutate, path) in cases
                err = try
                    configuration_from(withkey(mutate))
                    nothing
                catch e
                    e
                end
                @test err isa ConfigurationError
                @test err.path == path
            end
        end

        @testset "cross-key constraints" begin
            # The receptor spacing must fit inside the grid.
            err = try
                configuration_from(withkey(t -> (t["grid"]["spacing"] = t["grid"]["extent"])))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "grid.spacing"

            # The high deposition velocity cannot sit below the low one.
            err = try
                configuration_from(
                    withkey(t -> (t["nuclide"]["deposition_velocity_high"] = 1e-6)),
                )
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "nuclide.deposition_velocity_high"

            # The washout table is defined only at tabulated intensities.
            err = try
                configuration_from(withkey(t -> (t["precipitation"]["rate"] = 2.0)))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "precipitation.rate"
        end

        @testset "the wind rose is checked on the way in" begin
            # One frequency per sector.
            err = try
                configuration_from(
                    withkey(t -> (t["wind_rose"]["frequencies"] = fill(0.125, 8))),
                )
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "wind_rose.frequencies"

            # One stability fraction per Pasquill class.
            err = try
                configuration_from(withkey(t -> (t["wind_rose"]["stability"] = fill(0.25, 4))))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "wind_rose.stability"

            # Sector count must be even and at least four.
            for n in (3, 15, 2)
                err = try
                    configuration_from(withkey(t -> (t["wind_rose"]["sectors"] = n)))
                    nothing
                catch e
                    e
                end
                @test err isa ConfigurationError
                @test err.path == "wind_rose.sectors"
            end

            # Raw counts rather than a distribution reach the WindRose constructor.
            @test_throws ArgumentError configuration_from(
                withkey(t -> (t["wind_rose"]["frequencies"] = fill(10.0, 16))),
            )
        end

        @testset "the declared convention reaches the field" begin
            from = configuration_from(
                withkey(t -> (t["wind_rose"]["convention"] = "blowing_from")),
            )
            toward = configuration_from(
                withkey(t -> (t["wind_rose"]["convention"] = "blowing_toward")),
            )
            g = grid(from.rose)
            for k = 1:16
                @test frequency_toward(from.rose, k) ==
                      frequency_toward(toward.rose, opposite(g, k))
            end
            χ(c, k) = dilution_long_term(
                5000sin(sector_bearing(g, k)),
                5000cos(sector_bearing(g, k)),
                c.site,
                c.rose,
            )
            @test argmax([χ(from, k) for k = 1:16]) == opposite(g, argmax([χ(toward, k) for k = 1:16]))
        end

        @testset "types are enforced" begin
            err = try
                configuration_from(withkey(t -> (t["source"]["height"] = "tall")))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "source.height"
            @test occursin("expected", err.message)
            # An unadorned TOML integer is accepted where a float is wanted.
            @test configuration_from(withkey(t -> (t["source"]["height"] = 50))).site.source.height ==
                  50.0
        end
    end

    # Physics validation: the solver is checked against the analytic properties
    # of the Gaussian plume rather than only against itself. These are the
    # statements that would fail first if a normalisation, a reflection term or
    # a sector width were wrong.
    @testset "Physics validation" begin
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
        distances = (500.0, 2000.0, 10_000.0)

        # Integrating the three-dimensional field across the wind leaves the
        # crosswind-integrated concentration, 2 exp(−H²/2Σ_z²)/(√(2π) Σ_z u).
        # This tests the 2πΣ_yΣ_z normalisation and the ground reflection at once.
        @testset "crosswind integral" begin
            for class in PASQUILL_CLASSES, r in distances
                f(y) = dilution_instantaneous(y, -r, 0.0, site, class, 0.0)
                numeric, _ = quadgk(f, -3e4, 3e4; rtol = 1e-10)
                H = effective_height(r, site, class)
                Σz = corrected_vertical_dispersion(r, site, class)
                u = transport_wind_speed(site, class)
                analytic = 2 * exp(-H^2 / (2Σz^2)) / (sqrt(2π) * Σz * u)
                @test numeric ≈ analytic rtol = 1e-8
            end
        end

        # Nothing is created or destroyed: the released activity crossing any
        # downwind plane, carried at the transport speed, is the whole release.
        @testset "mass conservation" begin
            for class in PASQUILL_CLASSES, r in distances
                u = transport_wind_speed(site, class)
                inner(z) = quadgk(
                    y -> dilution_instantaneous(y, -r, z, site, class, 0.0),
                    -3e4,
                    3e4;
                    rtol = 1e-9,
                )[1]
                total, _ = quadgk(inner, 0.0, 5000.0; rtol = 1e-8)
                @test u * total ≈ 1 rtol = 1e-6
            end
        end

        # The sector-averaged form is the crosswind integral spread uniformly
        # over the arc the sector subtends, since 2/√(2π) = √(2/π). Any error in
        # the sector width or in that constant shows up here.
        @testset "sector average is the crosswind integral over the arc" begin
            for sectors in (SectorGrid(8), SectorGrid(16), SectorGrid(36))
                θ_L = sector_width(sectors)
                for class in PASQUILL_CLASSES, r in distances
                    extended = dilution_extended(0.0, -r, site, class, 0.0, sectors)
                    H = effective_height(r, site, class)
                    Σz = corrected_vertical_dispersion(r, site, class)
                    u = transport_wind_speed(site, class)
                    crosswind = 2 * exp(-H^2 / (2Σz^2)) / (sqrt(2π) * Σz * u)
                    @test extended ≈ crosswind / (r * θ_L) rtol = 1e-10
                end
            end
        end

        # The ground-level centreline maximum falls at Σ_z = H/√2. That is exact
        # only under the assumptions of its derivation — H fixed and Σ_y ∝ Σ_z —
        # so it is asserted there. In the real field it lands within 6 %, and the
        # residual is Σ_y/Σ_z drifting by about a quarter across the peak region,
        # not an error.
        @testset "ground-level maximum" begin
            for H in (50.0, 100.0, 200.0), u in (2.0, 5.0)
                f(σz) = exp(-H^2 / (2σz^2)) / (π * (2σz) * σz * u)
                σzs = range(1.0, 5H, length = 200_000)
                @test σzs[argmax(f.(σzs))] ≈ H / sqrt(2) rtol = 1e-4
            end
            xs = 10 .^ range(1.5, 5, length = 4000)
            χ = [dilution_instantaneous(0.0, -x, 0.0, site, PASQUILL_D, 0.0) for x in xs]
            xm = xs[argmax(χ)]
            ratio =
                corrected_vertical_dispersion(xm, site, PASQUILL_D) /
                (effective_height(xm, site, PASQUILL_D) / sqrt(2))
            @test 0.9 < ratio < 1.1
        end
    end

    # Validation against the published literature, as distinct from the analytic
    # self-consistency above. Several of the parameterisations this code inherits
    # from the normative turn out to be standard published schemes, and where
    # they are, the identity is asserted rather than described.
    @testset "Literature validation" begin
        # The vertical dispersion scheme is Hosker's fit to F.B. Smith (1972)
        # and Briggs (1973), published as IAEA-SM-181/19 (1974) and printed in
        # Smith and Simmonds (eds.), HPA-RPD-058, Health Protection Agency
        # (2009), Table 3.3, and in Clarke, NRPB-R91 (1979), Table 3. Both
        # printings agree digit for digit, and the ORNL codes that implement it
        # (ORNL-5913 Table 6, ORNL/TM-6874 p. 25) carry the same numbers.
        #
        # Transcribing a table is exactly where a package of tabulated constants
        # fails, so the whole table is asserted rather than sampled.
        @testset "vertical dispersion coefficients are Hosker's" begin
            # σ_z = a x^b / (1 + c x^d), x and σ_z in metres.
            shape = (
                (0.112, 1.06, 5.38e-4, 0.815),    # A
                (0.130, 0.950, 6.52e-4, 0.750),   # B
                (0.112, 0.920, 9.05e-4, 0.718),   # C
                (0.098, 0.889, 1.35e-3, 0.688),   # D
                (0.0609, 0.895, 1.96e-3, 0.684),  # E
                (0.0638, 0.783, 1.36e-3, 0.672),  # F
            )
            for (i, class) in enumerate(PASQUILL_CLASSES)
                c = vertical_shape_coefficients(class)
                @test (c.a₁, c.b₁, c.a₂, c.b₂) == shape[i]
            end

            # F(z₀,x) = ln(f x^g [1 + {h x^j}⁻¹]) for z₀ > 0.1 m,
            # F(z₀,x) = ln(f x^g [1 + h x^j]⁻¹)   for z₀ ≤ 0.1 m.
            correction = (
                (0.01, 1.56, 0.0480, 6.25e-4, 0.45),
                (0.04, 2.02, 0.0269, 7.76e-4, 0.37),
                (0.10, 2.72, 0.0, 0.0, 0.0),
                (0.40, 5.16, -0.098, 18.6, -0.225),
                (1.00, 7.37, -0.0957, 4.29e3, -0.60),
                (4.00, 11.7, -0.128, 4.59e4, -0.78),
            )
            for (i, roughness) in enumerate(ROUGHNESS_CLASSES)
                c = roughness_coefficients(roughness)
                @test (c.z₀, c.c₁, c.d₁, c.c₂, c.d₂) == correction[i]
            end

            # z₀ = 10 cm is the reference the fit is normalised on, so its
            # correction is ln(2.72) — unity to within the rounding of e.
            @test roughness_correction(1000.0, ROUGHNESS_PASTURE) ≈ 1 atol = 7e-4
        end

        # The sector-averaged long-term form is a regulatory equation, stated
        # identically by NRC Regulatory Guide 1.111 Rev. 1 (1977) Eq. (3) and
        # XOQDOQ, NUREG/CR-2919 (1982) Eq. (1), by IAEA Safety Reports Series
        # No. 19 (2001) Eq. (V-2), and by the German AVV zu §47 StrlSchV (2012)
        # Eq. (4.4). RG 1.111 writes the constant as 2.032 and states in words
        # that it is √(2/π) divided by a 22.5° sector in radians; SRS-19 works
        # in twelve sectors, where the same constant is 1.5238.
        @testset "the sector constant is the published regulatory value" begin
            # SRS-19 Eq. (3) prints the twelve-sector constant as 12/√(2π³),
            # which is a closed form rather than a rounded decimal; the
            # sixteen-sector analogue is RG 1.111's 2.032.
            @test sqrt(2 / π) * 16 / (2π) ≈ 16 / sqrt(2π^3) rtol = 1e-14
            @test sqrt(2 / π) * 12 / (2π) ≈ 12 / sqrt(2π^3) rtol = 1e-14
            @test sqrt(2 / π) * 16 / (2π) ≈ 2.032 rtol = 2e-4
            @test sqrt(2 / π) * 12 / (2π) ≈ 1.5238 rtol = 2e-4

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
            distances = (500.0, 2000.0, 10_000.0)

            for (n, published) in ((16, 2.032), (12, 1.5238))
                sectors = SectorGrid(n)
                for class in PASQUILL_CLASSES, r in distances
                    H = effective_height(r, site, class)
                    Σz = corrected_vertical_dispersion(r, site, class)
                    u = transport_wind_speed(site, class)
                    regulatory = published * exp(-H^2 / (2Σz^2)) / (Σz * u * r)
                    @test dilution_extended(0.0, -r, site, class, 0.0, sectors) ≈ regulatory rtol =
                        2e-4
                end
            end
        end

        # IAEA Safety Series No. 57, Generic Models and Parameters for Assessing
        # the Environmental Transfer of Radionuclides from Routine Releases,
        # STI/PUB/611, Vienna (1982), §3.6, Eq. (3.14A). Read in the primary
        # document, which states verbatim: "Reference values for the terms in
        # Eq. (3.14A) are 10⁻⁵ m⁻¹, 10⁻⁹ m⁻¹ for A and B respectively and
        # 1 × 10⁻² d⁻¹ and 2 × 10⁻⁵ d⁻¹ for λ₁ and λ₂ respectively."
        #
        # Safety Series 57 has been superseded — every page of it now carries a
        # "no longer valid" stamp — but it is the source these constants came
        # from, by way of reference [3] of CNCAN NSR-23, and no successor
        # restates them. See RESUSPENSION_MAXWELL_ANSPAUGH for the modern form.
        @testset "resuspension is IAEA Safety Series 57" begin
            c = RESUSPENSION_COEFFICIENTS
            @test c.A == 1e-5
            @test c.B == 1e-9
            @test c.λ₁ == 1e-2
            @test c.λ₂ == 2e-5

            # The same paragraph brackets them across the literature.
            @test 1e-6 ≤ c.A ≤ 1e-4
            @test 1e-10 ≤ c.B ≤ 1e-8
            # "of the order of weeks" and "in the range 50 to 100 years"
            @test 14 ≤ log(2) / c.λ₁ ≤ 120
            @test 50 ≤ log(2) / c.λ₂ / 365.25 ≤ 100

            @test resuspension_factor(0) ≈ c.A + c.B
            @test resuspension_factor(1e9) ≈ 0 atol = 1e-12
        end

        # The washout table follows the published intensity dependence and
        # brackets the published amplitudes; it does not reproduce a single
        # tabulation, and is asserted as what it is.
        @testset "washout follows the published intensity law" begin
            # Λ ∝ J^0.75, from Slinn (1977) via NRPB-R322 (ADMLC, 2001) §3.1.1:
            # "a net dependence of scavenging coefficient on J^0.75".
            for precipitation in (PRECIPITATION_RAIN, PRECIPITATION_SNOW),
                field in (:low, :high)

                Λ = [
                    getfield(washout_coefficients(precipitation, j), field) for
                    j in PRECIPITATION_RATES
                ]
                x = log.(collect(PRECIPITATION_RATES))
                y = log.(Λ)
                x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
                p = sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
                @test 0.65 ≤ p ≤ 0.80
            end

            # At 1 mm/h the German AVV zu §47 StrlSchV (2012), Anhang 7
            # Tabelle 3, gives 7e-5 s⁻¹ for aerosols and elemental iodine and
            # 3.5e-5 for tritiated water. Both lie inside this table's range.
            w = washout_coefficients(PRECIPITATION_RAIN, 1.0)
            @test w.low ≤ 7.0e-5 ≤ w.high
            @test w.low ≤ 3.5e-5 ≤ w.high

            # IAEA Safety Series 57 Table II gives the washout coefficient as
            # Λ = a·I with a = 1.6e-4 h(mm·s)⁻¹ for particulates and 1.1e-4 for
            # elemental iodine, so 1.6e-4 and 1.1e-4 s⁻¹ at 1 mm/h. Both fall
            # inside both species rows of the normative's rain bracket.
            for species in (WASHOUT_TRITIUM_IODINE, WASHOUT_OTHER_NUCLIDES)
                r = washout_coefficients(
                    PRECIPITATION_RAIN,
                    1.0,
                    WASHOUT_NORMATIVE,
                    species,
                )
                @test r.low ≤ 1.6e-4 ≤ r.high
                @test r.low ≤ 1.1e-4 ≤ r.high
            end

            # Safety Series 57 makes the dependence linear in intensity, where
            # the normative's table fits 0.75 and NRPB-R322 publishes 0.75.
            # NRPB-R157 §D3.4 brackets the exponent at 0.5 to 1.0, which
            # contains both positions, so neither is asserted against the other.
            @test 0.5 ≤ 0.75 ≤ 1.0
            @test 0.5 ≤ 1.0 ≤ 1.0

            # The snow columns are the rain columns scaled down by a constant —
            # a particle-scavenging suppression, and not a measurement.
            for (field, factor) in ((:low, 100), (:high, 500))
                ratios = [
                    getfield(washout_coefficients(PRECIPITATION_RAIN, j), field) /
                    getfield(washout_coefficients(PRECIPITATION_SNOW, j), field) for
                    j in PRECIPITATION_RATES
                ]
                @test all(r -> isapprox(r, factor; rtol = 0.2), ratios)
            end
        end

        # The building-wake coefficient. IAEA Safety Reports Series No. 19
        # (2001) Eq. (6) writes the wake-broadened vertical parameter as
        # (σ_z² + A_B/π)^(1/2), and the German AVV zu §47 StrlSchV (2012)
        # Eqs. (4.31)/(4.32) as √(σ² + I_G²/π): a coefficient of one in both.
        @testset "wake coefficient is the published one" begin
            @test DEFAULT_WAKE_COEFFICIENT == 1.0
            @test NORMATIVE_WAKE_COEFFICIENT == 1.5

            b = Building(; east = 25.0, north = 0.0, height = 60.0, frontal_area = 3600.0)
            iaea = BuildingEnvelope([b])
            normative = BuildingEnvelope([b]; wake_coefficient = NORMATIVE_WAKE_COEFFICIENT)
            A = equivalent_area(iaea)
            for σ in (10.0, 50.0, 200.0)
                # released inside the cavity, where the correction is undiluted
                @test wake_broadened(σ, 0.0, iaea) ≈ sqrt(σ^2 + A / π)
                @test wake_broadened(σ, 0.0, normative) ≈ sqrt(σ^2 + 1.5A / π)
                @test wake_broadened(σ, 0.0, iaea) < wake_broadened(σ, 0.0, normative)
            end
        end

        # Ogram, Precipitation Scavenging of Tritiated Water Vapour (HTO),
        # Ontario Hydro Research Division 85-233-K (1985), §6.0 Eq. (38), with
        # its Table V at 0.5, 1 and 2 mm/h. Snow scavenging of HTO is isotopic
        # exchange at the crystal surface, so the particle suppression the
        # normative applies is the wrong physics for it.
        @testset "HTO snow washout is Ogram 1985" begin
            @test ogram_snow_washout(0.5) ≈ 2.9e-4 rtol = 0.02
            @test ogram_snow_washout(1.0) ≈ 4.2e-4 rtol = 0.02
            @test ogram_snow_washout(2.0) ≈ 6.1e-4 rtol = 0.02
            @test ogram_snow_washout(1.0) ≈ 1.2e-4 + 3.0e-4

            # Ogram's measurement lands on NSR-23's own other-nuclides snow
            # lower limit, a row it does not cite Ogram for — which is the
            # evidence that the tritium row's snow column is the anomaly.
            for r in PRECIPITATION_RATES
                other = washout_coefficients(
                    PRECIPITATION_SNOW,
                    r,
                    WASHOUT_NORMATIVE,
                    WASHOUT_OTHER_NUCLIDES,
                )
                @test 0.8 < ogram_snow_washout(r) / other.low < 1.2
            end

            # NSR-23 Table 7 has two species rows; the other-nuclides snow
            # values run three to four orders of magnitude above the tritium
            # row's, in the opposite direction to the particle suppression.
            for r in PRECIPITATION_RATES
                trit = washout_coefficients(
                    PRECIPITATION_SNOW,
                    r,
                    WASHOUT_NORMATIVE,
                    WASHOUT_TRITIUM_IODINE,
                )
                other = washout_coefficients(
                    PRECIPITATION_SNOW,
                    r,
                    WASHOUT_NORMATIVE,
                    WASHOUT_OTHER_NUCLIDES,
                )
                @test other.low / trit.low > 1e3
            end

            # rain is untouched by the choice; snow differs by about 10³
            for r in PRECIPITATION_RATES
                @test washout_coefficients(PRECIPITATION_RAIN, r, WASHOUT_HTO) ==
                      washout_coefficients(PRECIPITATION_RAIN, r, WASHOUT_NORMATIVE)
                ratio =
                    washout_coefficients(PRECIPITATION_SNOW, r, WASHOUT_HTO).high /
                    washout_coefficients(PRECIPITATION_SNOW, r, WASHOUT_NORMATIVE).high
                @test 900 < ratio < 1500
            end
        end

        # The combined momentum-and-buoyancy law must reduce to the pure
        # buoyancy law when the momentum flux vanishes. Briggs writes its
        # buoyancy term as 3F x²/(2β²u³) with β = 0.6, a denominator of
        # 2β² = 0.72, and that is now the default; the 2021 code had 0.5, which
        # overshoots the two-thirds law by 13.6 % and is kept as THESIS_RISE.
        @testset "combined rise reduces to the two-thirds law" begin
            # With c = 2β² = 0.72 the combined law implies a two-thirds
            # coefficient of (3/0.72)^(1/3) = 1.60915, against the 1.6 the
            # literature rounds it to and that `buoyant_rise` uses. The two
            # therefore agree to 0.6 %, not exactly — the same rounding the
            # neutral final rise shows as 21.425 against a published 21.4.
            @test (3 / BRIGGS_RISE.combined_buoyancy)^(1 / 3) ≈ 1.60915 rtol = 1e-4

            checked = 0
            for F in (5.0, 50.0, 500.0), u in (2.0, 5.0, 9.0), x in (50.0, 100.0, 300.0)
                two_thirds = 1.6 * F^(1 / 3) * x^(2 / 3) / u
                # The combined law is capped at the sum of the two final rises,
                # and the cap differs between the two coefficient sets, so only
                # assert where neither is binding.
                cap(r) =
                    final_momentum_rise(1e-14, 10.0, 2.0, u, -1e-6, r) +
                    final_buoyant_rise(F, u, -1e-6, r)
                briggs = combined_rise(x, F, 1e-14, 10.0, u, -1e-6, 2.0, BRIGGS_RISE)
                thesis = combined_rise(x, F, 1e-14, 10.0, u, -1e-6, 2.0, THESIS_RISE)
                (briggs < 0.99cap(BRIGGS_RISE) && thesis < 0.99cap(THESIS_RISE)) || continue
                checked += 1
                @test briggs ≈ two_thirds rtol = 6e-3
                @test thesis / two_thirds ≈ (6 / 1.6^3)^(1 / 3) rtol = 2e-3
            end
            @test checked ≥ 6          # the sweep must actually exercise the branch
            @test BRIGGS_RISE.combined_buoyancy == 2 * 0.6^2       # Briggs, β = 0.6
            @test 3 / 1.6^3 ≈ 0.7324 rtol = 1e-4                   # the same from the law
            @test THESIS_RISE.combined_buoyancy == 0.5
        end

        # The neutral momentum rise is 3 w₀D/u in Briggs (1969) Eq. 5.2, as
        # implemented by EPA ISC3 Eq. (1-16). The 2021 code had 1.5, unsourced.
        @testset "neutral momentum rise is Briggs" begin
            for w₀ in (5.0, 15.0), D in (1.0, 3.0), u in (2.0, 8.0)
                @test final_momentum_rise(0.0, w₀, D, u, -1e-6, BRIGGS_RISE) ≈
                      3 * w₀ * D / u
                @test final_momentum_rise(0.0, w₀, D, u, -1e-6, THESIS_RISE) ≈
                      1.5 * w₀ * D / u
            end
            @test BRIGGS_RISE.neutral_momentum == 3.0
        end

        # The stable final rise coefficient is 2.6 in Briggs and the Handbook on
        # Atmospheric Diffusion, and 2.4 in NRC XOQDOQ.
        @testset "stable rise coefficient, Briggs against XOQDOQ" begin
            @test BRIGGS_RISE.stable_final == 2.6
            @test XOQDOQ_RISE.stable_final == 2.4
            for F in (5.0, 200.0), u in (3.0, 8.0), S in (1e-4, 1e-3)
                b = final_buoyant_rise(F, u, S, BRIGGS_RISE)
                x = final_buoyant_rise(F, u, S, XOQDOQ_RISE)
                @test x ≤ b
            end
        end

        # Briggs (1973) open-country lateral dispersion, as tabulated in Hanna,
        # Briggs and Hosker, Handbook on Atmospheric Diffusion, DOE/TIC-11223
        # (1982), Table 4.5, attributing ATDL Contribution No. 79:
        # σ_y = a x (1 + 10⁻⁴x)^(−1/2) with a running
        # 0.22, 0.16, 0.11, 0.08, 0.06, 0.04 from class A to F. The table states
        # its own validity band as 10² < x < 10⁴ m.
        @testset "σ_y is Briggs 1973 open country" begin
            briggs_a = (0.22, 0.16, 0.11, 0.08, 0.06, 0.04)
            for (i, class) in enumerate(PASQUILL_CLASSES)
                @test lateral_coefficient(class) == briggs_a[i]
                for x in (10.0, 100.0, 1000.0, 10_000.0, 50_000.0)
                    @test lateral_dispersion(x, class) ≈
                          briggs_a[i] * x * (1 + 1e-4 * x)^(-0.5) rtol = 1e-12
                end
            end
        end

        # Briggs distance to final rise: x_f = 14F^(5/8) below 55 m⁴/s³,
        # 34F^(2/5) above.
        @testset "distance to final rise is Briggs" begin
            for F in (1.0, 10.0, 54.9)
                @test buoyancy_transition_distance(F) ≈ 14 * F^(5 / 8) rtol = 1e-12
            end
            for F in (55.0, 200.0, 5000.0)
                @test buoyancy_transition_distance(F) ≈ 34 * F^(2 / 5) rtol = 1e-12
            end
        end

        # The neutral final buoyant rise is written here as
        # 1.6F^(1/3)(3.5x_f)^(2/3)/u. Briggs publishes it as 21.4F^(3/4)/u and
        # 38.7F^(3/5)/u. Substituting x_f shows these are one expression:
        # 1.6·49^(2/3) = 21.425 and 1.6·119^(2/3) = 38.71, which the literature
        # rounds. The agreement is therefore exact up to that rounding, and the
        # test asserts both the algebra and the numbers.
        @testset "neutral final rise is the published Briggs form" begin
            @test 1.6 * 49^(2 / 3) ≈ 21.425 rtol = 1e-4
            @test 1.6 * 119^(2 / 3) ≈ 38.71 rtol = 1e-4
            for F in (5.0, 20.0, 54.0), u in (3.0, 8.0)
                @test final_buoyant_rise(F, u, -1e-6) ≈ 21.4 * F^0.75 / u rtol = 2e-3
            end
            for F in (56.0, 200.0, 1000.0), u in (3.0, 8.0)
                @test final_buoyant_rise(F, u, -1e-6) ≈ 38.7 * F^0.6 / u rtol = 2e-3
            end
        end

        # Briggs stable final rise, 2.6[F/(u s)]^(1/3), and the two-thirds law
        # for the transitional rise.
        @testset "stable final rise and the two-thirds law" begin
            for F in (10.0, 100.0), u in (2.0, 6.0), S in (1e-4, 1e-3)
                stable = 2.6 * (F / (u * S))^(1 / 3)
                @test final_buoyant_rise(F, u, S) <= stable * (1 + 1e-12)
                # where the stable branch is the binding one, it is exact
                if stable <
                   1.6 * F^(1/3) * (3.5 * buoyancy_transition_distance(F))^(2/3) / u &&
                   stable < 5.0 * F^(1/4) * S^(-3/8)
                    @test final_buoyant_rise(F, u, S) ≈ stable rtol = 1e-12
                end
            end
            for F in (10.0, 100.0), u in (2.0, 6.0), x in (10.0, 100.0)
                two_thirds = 1.6 * F^(1 / 3) * x^(2 / 3) / u
                r = buoyant_rise(x, F, u, 1e-3)
                @test r ≈ min(two_thirds, final_buoyant_rise(F, u, 1e-3)) rtol = 1e-12
            end
        end

        # σ_z is *not* Briggs — it is the g(x)F(x) form of the normative, an
        # NRPB-R91-style scheme with an explicit roughness correction. It is not
        # asserted equal to Briggs; it is asserted to agree within the factor
        # that separates published σ schemes, with the expected sign: less
        # vertical spread than Briggs in unstable air, more in stable.
        @testset "σ_z brackets Briggs open country" begin
            briggs_z = Dict(
                PASQUILL_A => (x -> 0.20x),
                PASQUILL_B => (x -> 0.12x),
                PASQUILL_C => (x -> 0.08x * (1 + 0.0002x)^(-0.5)),
                PASQUILL_D => (x -> 0.06x * (1 + 0.0015x)^(-0.5)),
                PASQUILL_E => (x -> 0.03x * (1 + 0.0003x)^(-1)),
                PASQUILL_F => (x -> 0.016x * (1 + 0.0003x)^(-1)),
            )
            for class in PASQUILL_CLASSES, x in (100.0, 300.0, 1000.0, 3000.0, 10_000.0)
                ratio =
                    vertical_dispersion(x, class, ROUGHNESS_PASTURE) / briggs_z[class](x)
                @test 0.4 < ratio < 1.6
            end
            # the sign of the difference, which is systematic
            @test vertical_dispersion(1000.0, PASQUILL_A, ROUGHNESS_PASTURE) <
                  briggs_z[PASQUILL_A](1000.0)
            @test vertical_dispersion(1000.0, PASQUILL_F, ROUGHNESS_PASTURE) >
                  briggs_z[PASQUILL_F](1000.0)
        end
    end

    @testset "WindRose" begin
        g = SectorGrid(16)
        uniform_stability = fill(1 / 6, 6)

        @testset "the convention is a half-turn rotation" begin
            f = zeros(16)
            f[1] = 1.0   # all the wind is in sector 1 (north)
            from = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
            toward = WindRose(g, f, BlowingToward(); stability = uniform_stability)

            # Wind from the north transports towards the south, sector 9.
            @test frequency_toward(from, 9) == 1.0
            @test frequency_toward(from, 1) == 0.0
            @test frequency_from(from, 1) == 1.0

            # Declared the other way, the same numbers mean the opposite field.
            @test frequency_toward(toward, 1) == 1.0
            @test frequency_toward(toward, 9) == 0.0
            @test frequency_from(toward, 9) == 1.0

            # The two readings of one table are exactly opposite everywhere.
            for k = 1:16
                @test frequency_toward(from, k) == frequency_toward(toward, opposite(g, k))
                @test frequency_from(from, k) == frequency_toward(from, opposite(g, k))
            end
        end

        @testset "a symmetric rose is invariant under the convention" begin
            f = fill(1 / 16, 16)
            a = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
            b = WindRose(g, f, BlowingToward(); stability = uniform_stability)
            @test all(k -> frequency_toward(a, k) == frequency_toward(b, k), 1:16)
        end

        @testset "stability fractions" begin
            f = fill(1 / 16, 16)
            s = [0.06533, 0.06533, 0.06533, 0.488, 0.158, 0.158]  # the 2021 input table
            rose = WindRose(g, f, BlowingFrom(); stability = s)
            for k = 1:16
                @test sum(stability_fraction(rose, k, c) for c in PASQUILL_CLASSES) ≈ 1 atol =
                    1e-3
                @test stability_fraction(rose, k, PASQUILL_D) == 0.488
            end

            # A per-sector matrix must follow the sector rotation of its rose.
            m = repeat(reshape(s, 1, 6), 16, 1)
            m[1, :] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            rose_m = WindRose(g, f, BlowingToward(); stability = m)
            @test stability_fraction(rose_m, 1, PASQUILL_D) == 1.0
            @test stability_fraction(rose_m, 2, PASQUILL_D) == 0.488
        end

        @testset "validation" begin
            f = fill(1 / 16, 16)
            @test_throws DimensionMismatch WindRose(
                g,
                fill(1 / 8, 8),
                BlowingFrom();
                stability = uniform_stability,
            )
            @test_throws DimensionMismatch WindRose(
                g,
                f,
                BlowingFrom();
                stability = fill(1 / 5, 5),
            )
            @test_throws DimensionMismatch WindRose(
                g,
                f,
                BlowingFrom();
                stability = fill(1 / 6, 8, 6),
            )
            # Raw counts rather than a normalised distribution.
            @test_throws ArgumentError WindRose(
                g,
                fill(10.0, 16),
                BlowingFrom();
                stability = uniform_stability,
            )
            bad = copy(f)
            bad[1] = -bad[2]
            bad[2] *= 2
            @test_throws ArgumentError WindRose(
                g,
                bad,
                BlowingFrom();
                stability = uniform_stability,
            )
            nan = copy(f)
            nan[1] = NaN
            @test_throws ArgumentError WindRose(
                g,
                nan,
                BlowingFrom();
                stability = uniform_stability,
            )
            @test_throws ArgumentError WindRose(
                g,
                f,
                BlowingFrom();
                stability = fill(1 / 3, 6),
            )
        end

        @testset "accessors" begin
            f = fill(1 / 16, 16)
            rose = WindRose(g, f, BlowingFrom(); stability = uniform_stability)
            @test grid(rose) === g
            @test nsectors(rose) == 16
            @test_throws BoundsError frequency_toward(rose, 0)
            @test_throws BoundsError stability_fraction(rose, 17, PASQUILL_A)
            @test occursin("WindRose", sprint(show, MIME"text/plain"(), rose))
        end
    end
end
