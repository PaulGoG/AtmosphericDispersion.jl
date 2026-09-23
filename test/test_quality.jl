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
# stay concretely typed, or the whole of QuadGK resolves dynamically. The
# integrand closes over the site, so both site types are guarded.
@testset "Inference at the quadrature boundary" begin
    site = REFERENCE_SITE
    plume = PrescribedPlume(site; height = 80.0, wind = 5.0)
    for s in (site, plume)
        @test Base.return_types(depletion_integral, (Float64, typeof(s), PasquillClass)) ==
              [Float64]
        @test @inferred(depletion_integral(1000.0, s, PASQUILL_D)) isa Float64
        @test @inferred(dry_depletion_factor(1000.0, s, PASQUILL_D, TRITIATED_WATER)) isa
              Float64
    end
    @test @inferred(effective_height(1000.0, site, PASQUILL_D)) isa Float64
    @test @inferred(corrected_vertical_dispersion(1000.0, site, PASQUILL_D)) isa Float64
end

# The kernels are written against AbstractSite. A PrescribedPlume holds its
# parameters as Union{Nothing,Float64}, and that union must not reach the return
# type of anything the kernels call.
@testset "Inference of the site interface" begin
    site = REFERENCE_SITE
    prescribed = PrescribedPlume(
        site; height = 80.0, wind = 5.0, lateral = 120.0, vertical = 60.0,)
    for s in (site, PrescribedPlume(site), prescribed)
        @test @inferred(effective_height(1000.0, s, PASQUILL_D)) isa Float64
        @test @inferred(transport_wind_speed(s, PASQUILL_D)) isa Float64
        @test @inferred(corrected_lateral_dispersion(1000.0, s, PASQUILL_D)) isa Float64
        @test @inferred(corrected_vertical_dispersion(1000.0, s, PASQUILL_D)) isa Float64
        @test @inferred(mixing_layer(s)) isa MixingLayer
    end
end
