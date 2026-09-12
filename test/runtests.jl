using AtmosphericDispersion
using Test
using Aqua

@testset "AtmosphericDispersion.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(AtmosphericDispersion)
    end
    # Write your tests here.
end
