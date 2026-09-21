#!/usr/bin/env julia
#
# Long-term dispersion field for a configured run.
#
#     julia scripts/run.jl config/reference.toml
#
# Writes nothing; prints the sector maxima and the most exposed sector.

include(joinpath(@__DIR__, "activate.jl"))

using AtmosphericDispersion
using Printf

function main(args)
    length(args) == 1 || error("usage: julia scripts/run.jl <config.toml>")
    config = load_configuration(only(args))
    site = config.site
    rose = config.rose
    g = grid(rose)

    @printf(
        "Release          %.3g Bq of %s over %.3g s\n",
        config.activity,
        config.nuclide.name,
        config.release_duration
    )
    @printf(
        "Release height   %.1f m, effective %.1f m at 1 km (class D)\n",
        release_height(site),
        effective_height(1000.0, site, PASQUILL_D)
    )
    println()
    @printf(
        "%-5s  %12s  %14s  %14s\n",
        "",
        "chi/Q [s/m3]",
        "chi [Bq s/m3]",
        "dry dep [Bq/m2]"
    )

    best, bestk = 0.0, 1
    for k = 1:nsectors(g)
        β = sector_bearing(g, k)
        r = config.extent
        east, north = r * sin(β), r * cos(β)
        χQ = dilution_long_term(
            east,
            north,
            site,
            rose;
            nuclide = config.nuclide,
            washout_duration = config.washout_duration,
            precipitation = config.precipitation,
            rate = config.precipitation_rate,
            washout_model = config.washout_model,
        )
        χ = χQ * config.activity
        @printf(
            "%-5s  %12.4g  %14.4g  %14.4g\n",
            sector_name(g, k),
            χQ,
            χ,
            dry_deposition(χ, config.nuclide)
        )
        if χQ > best
            best, bestk = χQ, k
        end
    end

    println()
    @printf(
        "Most exposed sector at %.0f m: %s (wind from %s)\n",
        config.extent,
        sector_name(g, bestk),
        sector_name(g, opposite(g, bestk))
    )
    return nothing
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
