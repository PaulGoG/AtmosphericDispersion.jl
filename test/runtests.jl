using AtmosphericDispersion
using Test
using Aqua
using ExplicitImports
using JET
using QuadGK
using TOML

include("fixtures.jl")

@testset "AtmosphericDispersion.jl" begin
    include("test_quality.jl")
    include("test_sectors.jl")
    include("test_stability.jl")
    include("test_windrose.jl")
    include("test_tables.jl")
    include("test_dispersion.jl")
    include("test_source.jl")
    include("test_plumerise.jl")
    include("test_buildings.jl")
    include("test_site.jl")
    include("test_mixing.jl")
    include("test_dilution.jl")
    include("test_depletion.jl")
    include("test_deposition.jl")
    include("test_config.jl")
    include("test_physics_validation.jl")
    include("test_literature.jl")
end
