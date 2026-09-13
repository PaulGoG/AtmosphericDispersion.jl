using AtmosphericDispersion
using Documenter

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
    format = Documenter.HTML(;
        canonical = "https://PaulGoG.github.io/AtmosphericDispersion.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API" => [
            "Sectors, stability and the wind rose" => "api/geometry.md",
            "Site and dispersion parameters" => "api/parameters.md",
            "Dilution, depletion and deposition" => "api/fields.md",
            "Configuration" => "api/configuration.md",
        ],
    ],
)

deploydocs(; repo = "github.com/PaulGoG/AtmosphericDispersion.jl", devbranch = "main")
