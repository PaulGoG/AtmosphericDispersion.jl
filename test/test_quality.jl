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
    @test Base.return_types(depletion_integral, (Float64, Site, PasquillClass)) == [Float64]
    @test @inferred(depletion_integral(1000.0, site, PASQUILL_D)) isa Float64
    @test @inferred(dry_depletion_factor(1000.0, site, PASQUILL_D, TRITIATED_WATER)) isa
          Float64
    @test @inferred(effective_height(1000.0, site, PASQUILL_D)) isa Float64
    @test @inferred(corrected_vertical_dispersion(1000.0, site, PASQUILL_D)) isa Float64
end
