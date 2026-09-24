#!/usr/bin/env julia
#
# Long-term dispersion field for a configured run.
#
#     julia scripts/run.jl config/reference.toml
#
# Evaluates the long-term χ/Q on the receptor grid the configuration
# describes: every sector of the wind rose at the radii from `grid.spacing`
# to `grid.extent` in steps of `grid.spacing`. Prints, per sector, the
# dilution, the air concentration and the dry deposition at the extent and
# the maximum over the radii; then the most exposed receptor and the air
# concentration resuspended from its deposition under the configured
# `model.resuspension`. Writes nothing.

include(joinpath(@__DIR__, "activate.jl"))

using AtmosphericDispersion
using Printf

"Radii of the receptor grid, `spacing` to `extent` in steps of `spacing`, the extent included."
function grid_radii(config)
    radii = collect(config.spacing:config.spacing:config.extent)
    last(radii) < config.extent && push!(radii, config.extent)
    return radii
end

function main(args)
    length(args) == 1 || error("usage: julia scripts/run.jl <config.toml>")
    config = load_configuration(only(args))
    site = config.site
    rose = config.rose
    g = grid(rose)
    radii = grid_radii(config)

    @printf("Release          %.3g Bq of %s over %.3g s\n",
        config.activity,
        config.nuclide.name,
        config.release_duration)
    @printf("Release height   %.1f m, effective %.1f m at 1 km (class D)\n",
        release_height(site),
        effective_height(1000.0, site, PASQUILL_D))
    @printf("Receptors        %d sectors, %d radii from %.0f to %.0f m\n",
        nsectors(g),
        length(radii),
        first(radii),
        last(radii))
    scheme = dispersion_scheme(site)
    lower, upper = validity_range(scheme)
    outside = count(r -> !within_validity(r, scheme), radii)
    @printf("Dispersion       %s, fitted over %.0f to %.0f m; %d of %d radii extrapolate it\n",
        scheme,
        lower,
        upper,
        outside,
        length(radii))
    println()
    @printf("%-5s  %45s  %22s\n", "", "at the extent", "maximum over the radii")
    @printf("%-5s  %12s  %14s  %15s  %12s  %8s\n",
        "",
        "chi/Q [s/m3]",
        "chi [Bq s/m3]",
        "dry dep [Bq/m2]",
        "chi/Q [s/m3]",
        "at r [m]")

    best_extent = (0.0, 1)          # χ/Q at the extent, sector
    best_grid = (0.0, 1, 0.0)       # χ/Q, sector, radius
    for k in 1:nsectors(g)
        β = sector_bearing(g, k)
        profile = [dilution_long_term(
                       r * sin(β),
                       r * cos(β),
                       site,
                       rose;
                       nuclide = config.nuclide,
                       washout = config.washout,
                   ) for r in radii]
        χQ = last(profile)
        χ = χQ * config.activity
        i = argmax(profile)
        @printf("%-5s  %12.4g  %14.4g  %15.4g  %12.4g  %8.0f\n",
            sector_name(g, k),
            χQ,
            χ,
            dry_deposition(χ, config.nuclide),
            profile[i],
            radii[i])
        χQ > first(best_extent) && (best_extent = (χQ, k))
        profile[i] > first(best_grid) && (best_grid = (profile[i], k, radii[i]))
    end

    println()
    @printf("Most exposed sector at %.0f m: %s (wind from %s)\n",
        config.extent,
        sector_name(g, best_extent[2]),
        sector_name(g, opposite(g, best_extent[2])))
    χQ, k, r = best_grid
    deposition = dry_deposition(χQ * config.activity, config.nuclide)
    @printf("Most exposed receptor: %s at %.0f m, chi/Q %.4g s/m3, dry deposition %.4g Bq/m2\n",
        sector_name(g, k),
        r,
        χQ,
        deposition)
    @printf("Mean air concentration there over the release: %.4g Bq/m3\n",
        χQ * config.activity / config.release_duration)
    @printf("Resuspended from its deposition: %.4g Bq/m3 after 1 d, %.4g Bq/m3 after 1 a\n",
        resuspended_concentration(deposition, 1.0, config.resuspension),
        resuspended_concentration(deposition, 365.25, config.resuspension))
    return nothing
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main(ARGS)
end
