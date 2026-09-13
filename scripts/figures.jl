#!/usr/bin/env julia
#
# Representative figures for the README, from the reference configuration.
#
#     julia --project=scripts scripts/figures.jl
#
# Plotting lives in its own environment so that the package itself does not
# depend on Makie: a library that computes dispersion factors should not oblige
# every user of it to build a plotting stack.

using AtmosphericDispersion
using CairoMakie, MathTeXEngine
using Printf, Statistics

set_theme!(
    Theme(
        fonts = (;
            regular = texfont(:text),
            bold = texfont(:bold),
            italic = texfont(:italic),
        ),
        fontsize = 20,
        figure_padding = 16,
        Axis = (
            xgridstyle = :dash,
            ygridstyle = :dash,
            xgridcolor = (:grey, 0.12),
            ygridcolor = (:grey, 0.12),
            xminorticksvisible = false,
            yminorticksvisible = false,
            xtickalign = 1,
            ytickalign = 1,
        ),
    ),
)

"Okabe–Ito, colourblind-safe."
const PALETTE = (
    blue = "#0072B2",
    orange = "#E69F00",
    green = "#009E73",
    red = "#D55E00",
    purple = "#CC79A7",
    sky = "#56B4E9",
    black = "#000000",
)

const FIGURES = joinpath(@__DIR__, "..", "figures")

function reference_case()
    config = load_configuration(joinpath(@__DIR__, "..", "config", "reference.toml"))
    return config
end

"""
    dispersion_field(config; half_width, n)

Long-term dilution factor over a square receptor grid centred on the stack.

The sixteen-fold structure is the wind rose: the field is not circular, because
the frequency of transport towards each sector is not uniform. This is the
figure in which the direction convention is visible — reading the rose the
other way round would reflect the map through the origin.
"""
function dispersion_field(config; half_width = 15_000.0, n = 221)
    xs = range(-half_width, half_width, length = n)
    χ = Matrix{Float64}(undef, n, n)
    for (i, east) in enumerate(xs), (j, north) in enumerate(xs)
        χ[i, j] = dilution_long_term(
            east,
            north,
            config.site,
            config.rose;
            nuclide = config.nuclide,
        )
    end
    return xs, χ
end

function figure_field(config)
    xs, χ = dispersion_field(config)
    activity = config.activity
    # time-integrated concentration for the configured release, Bq s m^-3
    field = χ .* activity

    # Colour limits set from the field beyond 1 km. Inside that the elevated
    # plume has not reached the ground and the factor collapses through thirty
    # decades, which would take the whole colour scale and leave the sector
    # structure — the thing worth seeing — as one flat tone.
    far = [
        field[i, j] for (i, e) in enumerate(xs), (j, n) in enumerate(xs) if
        hypot(e, n) > 1000 && field[i, j] > 0
    ]
    lo, hi = quantile(far, 0.02), maximum(far)

    fig = Figure(size = (900, 720))
    ax = Axis(fig[1, 1], xlabel = "East [km]", ylabel = "North [km]", aspect = DataAspect())
    hm = heatmap!(
        ax,
        xs ./ 1000,
        xs ./ 1000,
        clamp.(field, lo, hi),
        colormap = :viridis,
        colorscale = log10,
        colorrange = (lo, hi),
    )
    scatter!(
        ax,
        [0.0],
        [0.0],
        color = :white,
        strokecolor = :black,
        strokewidth = 1.5,
        markersize = 13,
    )
    text!(ax, 0.4, 0.4; text = "Stack", fontsize = 15, color = :white)
    # Explicit ticks: the range spans about 1.2 decades, so automatic log ticks
    # land on fractional exponents like 10^6.75, which are not a thing anyone
    # reads off a colour bar.
    Colorbar(
        fig[1, 2],
        hm,
        label = L"Time-integrated concentration [Bq s m$^{-3}$]",
        ticks = (
            [2e6, 5e6, 1e7, 2e7],
            [L"2\times10^{6}", L"5\times10^{6}", L"10^{7}", L"2\times10^{7}"],
        ),
    )

    g = grid(config.rose)
    k = argmax([frequency_toward(config.rose, i) for i = 1:nsectors(g)])
    text!(
        ax,
        0.03,
        0.03;
        text = "Most exposed: " * sector_name(g, k),
        space = :relative,
        align = (:left, :bottom),
        fontsize = 16,
        color = :white,
    )

    path = joinpath(FIGURES, "dispersion_field.png")
    mkpath(FIGURES)
    save(path, fig; px_per_unit = 3)
    return path, sector_name(g, k)
end

function figure_regimes(config)
    site = config.site
    rose = config.rose
    class = PASQUILL_D
    r = 10 .^ range(log10(100.0), log10(30_000.0), length = 400)

    # wind from the north, so the plume runs south; sample along that axis
    inst = [dilution_instantaneous(0.0, -x, 0.0, site, class, 0.0) for x in r]
    ext = [dilution_extended(0.0, -x, site, class, 0.0) for x in r]
    lt = [dilution_long_term(0.0, -x, site, rose) for x in r]

    fig = Figure(size = (1150, 470))

    ax1 = Axis(
        fig[2, 1],
        xlabel = "Downwind distance [km]",
        ylabel = L"$\chi/Q$ [s m$^{-3}$]",
        xscale = log10,
        yscale = log10,
        xticks = ([0.1, 1.0, 10.0], ["0.1", "1", "10"]),
    )
    l1 = lines!(ax1, r ./ 1000, inst, color = PALETTE.blue, linewidth = 2)
    l2 = lines!(ax1, r ./ 1000, ext, color = PALETTE.orange, linewidth = 2)
    l3 = lines!(ax1, r ./ 1000, lt, color = PALETTE.green, linewidth = 2)
    # Clipped to the range that carries information. Below about 300 m the
    # elevated plume has not reached the ground and the factor falls through
    # thirty decades, which is physical and not worth thirty decades of axis.
    ylims!(ax1, 1e-12, 3e-7)

    ax2 = Axis(
        fig[2, 2],
        xlabel = "Downwind distance [km]",
        ylabel = "Effective release height [m]",
        xscale = log10,
        xticks = ([0.1, 1.0, 10.0], ["0.1", "1", "10"]),
    )
    handles = []
    colours = (
        PALETTE.blue,
        PALETTE.sky,
        PALETTE.green,
        PALETTE.orange,
        PALETTE.purple,
        PALETTE.red,
    )
    for (c, col) in zip(PASQUILL_CLASSES, colours)
        push!(
            handles,
            lines!(
                ax2,
                r ./ 1000,
                [effective_height(x, site, c) for x in r],
                color = col,
                linewidth = 1.8,
            ),
        )
    end
    hlines!(
        ax2,
        [site.source.height],
        color = PALETTE.black,
        linestyle = :dash,
        linewidth = 1.2,
    )
    text!(
        ax2,
        0.04,
        site.source.height + 2;
        text = "Stack height",
        fontsize = 14,
        align = (:left, :bottom),
    )

    Legend(
        fig[1, 1],
        [l1, l2, l3],
        ["Instantaneous, 3-D", "Extended, sector-averaged", "Long term, rose-weighted"],
        orientation = :horizontal,
        framevisible = false,
        labelsize = 14,
        nbanks = 2,
    )
    Legend(
        fig[1, 2],
        handles,
        ["A", "B", "C", "D", "E", "F"],
        orientation = :horizontal,
        framevisible = false,
        labelsize = 14,
        nbanks = 1,
        colgap = 10,
    )
    rowgap!(fig.layout, 1, 6)
    colgap!(fig.layout, 1, 30)

    path = joinpath(FIGURES, "regimes.png")
    mkpath(FIGURES)
    save(path, fig; px_per_unit = 3)
    return path
end

"""
    figure_animation(config)

The instantaneous ground-level plume, as the wind direction sweeps the compass.

This is the figure that shows the direction convention doing its work. The
bearing is the direction the wind blows **from**, as in a `.met` file, so the
plume always lies on the opposite side of the stack: at a bearing of 0 — a
northerly — the plume runs south. Averaging these over a year, weighted by how
often the wind blows towards each sector, is what the long-term field is.
"""
function figure_animation(config; half_width = 8_000.0, n = 141, frames = 72)
    xs = range(-half_width, half_width, length = n)
    bearings = range(0, 2π, length = frames + 1)[1:frames]
    site, class = config.site, PASQUILL_D

    field = Observable(zeros(n, n))
    label = Observable("")
    compute!(β) = begin
        f = [dilution_instantaneous(e, nn, 0.0, site, class, β) for e in xs, nn in xs]
        field[] = f
    end

    fig = Figure(size = (780, 640))
    ax = Axis(fig[1, 1], xlabel = "East [km]", ylabel = "North [km]", aspect = DataAspect())
    compute!(0.0)
    lo, hi = 1e-9, 3e-6
    hm = heatmap!(
        ax,
        xs ./ 1000,
        xs ./ 1000,
        lift(f -> clamp.(f, lo, hi), field),
        colormap = :viridis,
        colorscale = log10,
        colorrange = (lo, hi),
    )
    scatter!(
        ax,
        [0.0],
        [0.0],
        color = :white,
        strokecolor = :black,
        strokewidth = 1.5,
        markersize = 12,
    )
    text!(
        ax,
        0.03,
        0.96;
        text = label,
        space = :relative,
        align = (:left, :top),
        fontsize = 17,
        color = :white,
    )
    Colorbar(
        fig[1, 2],
        hm,
        label = L"$\chi/Q$ [s m$^{-3}$]",
        ticks = (
            [1e-9, 1e-8, 1e-7, 1e-6],
            [L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}"],
        ),
    )

    g = grid(config.rose)
    path = joinpath(FIGURES, "plume_sweep.gif")
    mkpath(FIGURES)
    record(fig, path, bearings; framerate = 12) do β
        compute!(β)
        label[] = "Wind from " * sector_name(g, sector_of(g, β))
    end
    return path
end

function main()
    config = reference_case()
    p1, sector = figure_field(config)
    println("wrote ", p1, "   (most exposed sector ", sector, ")")
    p2 = figure_regimes(config)
    println("wrote ", p2)
    p3 = figure_animation(config)
    println("wrote ", p3)
end

main()
