using AtmosphericDispersion
using Test
using Aqua
using ExplicitImports
using JET

@testset "AtmosphericDispersion.jl" begin
    @testset "Static quality assurance" begin
        Aqua.test_all(AtmosphericDispersion)
        @test check_no_implicit_imports(AtmosphericDispersion) === nothing
        @test check_no_stale_explicit_imports(AtmosphericDispersion) === nothing
        @test check_all_qualified_accesses_via_owners(AtmosphericDispersion) === nothing
        JET.test_package(AtmosphericDispersion)
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
            T2 = [  # Tabel_2.csv: z_0, c_1, d_1, c_2, d_2
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
                        @test vertical_dispersion(x, k, r) == original_σ_z(x, P, j)
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
