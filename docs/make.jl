include(joinpath(@__DIR__, "activate.jl"))

using AtmosphericDispersion
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

DocMeta.setdocmeta!(
    AtmosphericDispersion,
    :DocTestSetup,
    :(using AtmosphericDispersion);
    recursive = true,
)

makedocs(;
    modules = [AtmosphericDispersion],
    authors = "Paul-Adrian Gogîță",
    sitename = "AtmosphericDispersion.jl",
    plugins = [bib],
    format = Documenter.HTML(;
        canonical = "https://PaulGoG.github.io/AtmosphericDispersion.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Conventions" => "conventions.md",
        "Validation" => "validation.md",
        "The 2021 thesis code" => "thesis.md",
        "API" => [
            "Sectors, stability and the wind rose" => "api/geometry.md",
            "Site and dispersion parameters" => "api/parameters.md",
            "Dilution, depletion and deposition" => "api/fields.md",
            "Configuration" => "api/configuration.md",
        ],
        "References" => "references.md",
    ],
)

deploydocs(; repo = "github.com/PaulGoG/AtmosphericDispersion.jl", devbranch = "main")
