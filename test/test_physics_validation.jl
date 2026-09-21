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
    # Unbounded on purpose. The crosswind integral, the sector average and
    # the Σ_z = H/√2 maximum are all properties of the *unbounded* Gaussian;
    # under a lid the plume is confined and they no longer hold as written.
    # The lid has its own conservation test below.
    site = Site(; source = stack, atmosphere = air, mixing = MIXING_UNBOUNDED)
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
