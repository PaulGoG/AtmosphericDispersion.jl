#!/usr/bin/env julia
#
# Representative figures for the README, from the reference configuration.
#
#     julia scripts/figures.jl
#
# Plotting lives in its own environment so that the package itself does not
# depend on Makie: a library that computes dispersion factors should not oblige
# every user of it to build a plotting stack.

include(joinpath(@__DIR__, "activate.jl"))

using AtmosphericDispersion
using Printf, Statistics

include(joinpath(@__DIR__, "theme_common.jl"))

const FIGURES = joinpath(@__DIR__, "..", "figures")

# Both map axes carry the same title. Naming a direction in an axis title puts
# the word on the edge it sits against — "north" down the left-hand spine — and
# the compass marks state the orientation anyway.
const AXIS_LABEL = "Distance from the stack [km]"

"""
    compass!(ax; color = :white)

Mark the four cardinal directions on the edges of a map-frame axis.

An axis title naming a direction is read as marking the edge it sits on, which
puts north on the left and east along the bottom. The coordinate names go in the
axis titles and the directions go here, against the edge each one actually
points at.
"""
function compass!(ax; color = :white, fontsize = 19)
    marks = (
        ("N", 0.5, 0.985, (:center, :top)),
        ("S", 0.5, 0.015, (:center, :bottom)),
        ("E", 0.985, 0.5, (:right, :center)),
        ("W", 0.015, 0.5, (:left, :center)),
    )
    for (t, x, y, align) in marks
        text!(
            ax,
            x,
            y;
            text = t,
            space = :relative,
            align = align,
            fontsize = fontsize,
            font = :bold,
            color = color,
            # The marks sit over whatever the field happens to be doing at the
            # edge, which for a plume running due west is the bright end of the
            # colour scale.
            strokecolor = :black,
            strokewidth = 0.6,
        )
    end
    return ax
end

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
    ax = Axis(fig[1, 1], xlabel = AXIS_LABEL, ylabel = AXIS_LABEL, aspect = DataAspect())
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
    compass!(ax)
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
    # Clipped at the bottom only. Below about 300 m the elevated plume has not
    # reached the ground and the factor falls through thirty decades, which is
    # physical and not worth thirty decades of axis. The top must contain the
    # peaks: at 3e-7 the instantaneous and extended curves left the frame and
    # re-entered it, which reads as a break in the data.
    ylims!(ax1, 1e-12, 4e-6)

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
        0.11,
        site.source.height + 1.5;
        text = "Stack height",
        fontsize = 14,
        align = (:left, :bottom),
    )
    # Headroom below, so the stack-height line is not drawn on the frame.
    ylims!(ax2, site.source.height - 12, nothing)

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
function figure_animation(config; half_width = 20_000.0, n = 141, frames = 72)
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
    ax = Axis(fig[1, 1], xlabel = AXIS_LABEL, ylabel = AXIS_LABEL, aspect = DataAspect())
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
    compass!(ax)
    # Under the axes rather than over the field: the compass marks now occupy
    # the edges, and a white overlay at the top left ran into the N.
    Label(fig[2, 1], label, fontsize = 18, tellwidth = false)
    Colorbar(
        fig[1, 2],
        hm,
        label = L"$\chi/Q$ [s m$^{-3}$]",
        ticks = (
            [1e-9, 1e-8, 1e-7, 1e-6],
            [L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}"],
        ),
    )

    path = joinpath(FIGURES, "plume_sweep.gif")
    mkpath(FIGURES)
    record(fig, path, bearings; framerate = 12) do β
        compute!(β)
        # Both halves, deliberately. "Wind from N" alone is correct and reads
        # backwards at a glance — the plume is on the opposite side to the
        # direction named — and that ambiguity is what put the 2021 dose field
        # the wrong way round. The label states the consequence as well as the
        # convention so it cannot be misread.
        label[] = SWEEP_LABEL[](β)
    end
    return path
end

"""
    panel_sweep(path, panels; half_width, n, frames, columns, lo, hi, caption)

Record a multi-panel sweep of the compass, one panel per entry of `panels`.

Each panel is a `(title, f)` pair where `f(east, north, bearing)` returns the
ground-level dilution factor in s/m³. Every panel shares one colour scale, which
is the point: the comparison is quantitative, not a set of separately normalised
pictures.
"""
function panel_sweep(
    path,
    panels;
    half_width = 20_000.0,
    n = 101,
    frames = 48,
    columns = 3,
    lo = 1e-9,
    hi = 3e-6,
    caption = "",
    width = 1080,
    height = 780,
)
    xs = range(-half_width, half_width, length = n)
    bearings = range(0, 2π, length = frames + 1)[1:frames]
    fields = [Observable(zeros(n, n)) for _ in panels]
    label = Observable("")

    fig = Figure(size = (width, height))
    grid_layout = fig[1, 1] = GridLayout()
    rows = cld(length(panels), columns)
    # Ticks stop short of the frame. Labelling the extremes puts the last label
    # of one panel against the first of the next, which is what "2020" was.
    ticks = ([-10.0, 0.0, 10.0], ["−10", "0", "10"])
    axes = Axis[]
    local hm
    for (i, (title, f)) in enumerate(panels)
        row, col = fldmod1(i, columns)
        # The panel name goes above the frame. Inside it there is no corner
        # free: the compass occupies the midpoint of all four edges, and any
        # label long enough to be useful reaches one of them.
        ax = Axis(
            grid_layout[row, col],
            aspect = DataAspect(),
            title = title,
            titlesize = 16,
            xticks = ticks,
            yticks = ticks,
            xticklabelsvisible = row == rows,
            yticklabelsvisible = col == 1,
        )
        fields[i][] = [f(e, nn, 0.0) for e in xs, nn in xs]
        hm = heatmap!(
            ax,
            xs ./ 1000,
            xs ./ 1000,
            lift(v -> clamp.(v, lo, hi), fields[i]),
            colormap = :viridis,
            colorscale = log10,
            colorrange = (lo, hi),
        )
        # A smooth colour ramp over three decades hides everything but the
        # largest differences. Decade contours give the eye a fixed reference:
        # a panel where a contour sits further out is quantitatively different,
        # whatever the ramp looks like.
        contour!(
            ax,
            xs ./ 1000,
            xs ./ 1000,
            lift(v -> log10.(clamp.(v, lo, hi)), fields[i]),
            levels = collect(log10(lo):0.5:log10(hi)),
            color = (:white, 0.35),
            linewidth = 0.7,
        )
        scatter!(
            ax,
            [0.0],
            [0.0],
            color = :white,
            strokecolor = :black,
            strokewidth = 1.2,
            markersize = 8,
        )
        compass!(ax; fontsize = 13)
        # Top left, because the compass owns the top centre.
        push!(axes, ax)
    end

    # One label per direction for the whole grid rather than one per panel.
    Label(grid_layout[rows+1, 1:columns], AXIS_LABEL, fontsize = 17)
    Label(grid_layout[1:rows, 0], AXIS_LABEL, fontsize = 17, rotation = π / 2)

    Colorbar(
        fig[1, 2],
        hm,
        label = L"$\chi/Q$ [s m$^{-3}$]",
        ticks = (
            [1e-9, 1e-8, 1e-7, 1e-6],
            [L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}"],
        ),
    )
    Label(fig[2, 1:2], label, fontsize = 17, tellwidth = false)
    isempty(caption) || Label(fig[3, 1:2], caption, fontsize = 14, tellwidth = false)
    colgap!(grid_layout, 12)
    rowgap!(grid_layout, 12)
    rowgap!(fig.layout, 1, 4)
    length(fig.layout.content) > 3 && rowgap!(fig.layout, 2, 2)

    mkpath(FIGURES)
    record(fig, path, bearings; framerate = 12) do β
        for (i, (_, f)) in enumerate(panels)
            fields[i][] = [f(e, nn, β) for e in xs, nn in xs]
        end
        label[] = SWEEP_LABEL[](β)
    end
    return path
end

# Set once from the configured rose so every animation words the bearing the
# same way.
const SWEEP_LABEL = Ref{Function}(β -> "")

function set_sweep_label!(rose)
    g = grid(rose)
    SWEEP_LABEL[] =
        β ->
            "Wind from " *
            sector_name(g, sector_of(g, β)) *
            "  →  plume to " *
            sector_name(g, opposite(g, sector_of(g, β)))
    return nothing
end

"""
    comparison_sweep(path, base, variant; ...)

Two cases of the same release and the ratio between them, swept over the
compass.

The ratio panel is the point. Two dilution fields that differ by a factor of a
few look almost identical on a colour ramp spanning three decades, so the
difference is drawn explicitly on a diverging scale centred on 1, where anything
away from white is a real change and the colour bar reads as a factor.
"""
function comparison_sweep(
    path,
    base,
    variant;
    half_width = 12_000.0,
    n = 111,
    frames = 48,
    lo = 1e-9,
    hi = 3e-6,
    ratio_span = 10.0,
    caption = "",
    width = 1400,
    height = 500,
)
    xs = range(-half_width, half_width, length = n)
    bearings = range(0, 2π, length = frames + 1)[1:frames]
    fa, fb = Observable(zeros(n, n)), Observable(zeros(n, n))
    ratio = Observable(zeros(n, n))
    label = Observable("")

    compute!(β) = begin
        a = [base[2](e, nn, β) for e in xs, nn in xs]
        b = [variant[2](e, nn, β) for e in xs, nn in xs]
        fa[] = a
        fb[] = b
        # Only where both fields carry signal; elsewhere the ratio is noise
        # between two numbers that are both effectively zero.
        # eachindex over a matrix is linear, so a comprehension over it returns
        # a vector. Build in place to keep the shape.
        rv = similar(a)
        for i in eachindex(a, b)
            rv[i] =
                (a[i] > lo && b[i] > lo) ?
                clamp(b[i] / a[i], 1 / ratio_span, ratio_span) : 1.0
        end
        ratio[] = rv
    end
    compute!(0.0)

    fig = Figure(size = (width, height))
    grid_layout = fig[1, 1] = GridLayout()
    ticks = ([-10.0, 0.0, 10.0], ["−10", "0", "10"])
    local hm, hr
    for (i, (title, field)) in
        enumerate(((base[1], fa), (variant[1], fb), ("Ratio, with ÷ without", ratio)))
        ax = Axis(
            grid_layout[1, i],
            aspect = DataAspect(),
            title = title,
            titlesize = 16,
            xticks = ticks,
            yticks = ticks,
            yticklabelsvisible = i == 1,
        )
        if i < 3
            hm = heatmap!(
                ax,
                xs ./ 1000,
                xs ./ 1000,
                lift(v -> clamp.(v, lo, hi), i == 1 ? fa : fb),
                colormap = :viridis,
                colorscale = log10,
                colorrange = (lo, hi),
            )
            contour!(
                ax,
                xs ./ 1000,
                xs ./ 1000,
                lift(v -> log10.(clamp.(v, lo, hi)), i == 1 ? fa : fb),
                levels = collect(log10(lo):0.5:log10(hi)),
                color = (:white, 0.35),
                linewidth = 0.7,
            )
        else
            hr = heatmap!(
                ax,
                xs ./ 1000,
                xs ./ 1000,
                ratio,
                colormap = :RdBu,
                colorscale = log10,
                colorrange = (1 / ratio_span, ratio_span),
            )
        end
        scatter!(
            ax,
            [0.0],
            [0.0],
            color = :white,
            strokecolor = :black,
            strokewidth = 1.2,
            markersize = 8,
        )
        compass!(ax; fontsize = 13, color = i == 3 ? :black : :white)
    end
    Label(grid_layout[2, 1:3], AXIS_LABEL, fontsize = 16)
    Label(grid_layout[1, 0], AXIS_LABEL, fontsize = 16, rotation = π / 2)

    Colorbar(
        fig[1, 2],
        hm,
        label = L"$\chi/Q$ [s m$^{-3}$]",
        ticks = (
            [1e-9, 1e-8, 1e-7, 1e-6],
            [L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}"],
        ),
    )
    Colorbar(
        fig[1, 3],
        hr,
        label = "Ratio, second ÷ first",
        ticks = ([0.1, 0.3, 1.0, 3.0, 10.0], ["0.1", "0.3", "1", "3", "10"]),
    )
    Label(fig[2, 1:3], label, fontsize = 16, tellwidth = false)
    isempty(caption) || Label(fig[3, 1:3], caption, fontsize = 14, tellwidth = false)
    colgap!(grid_layout, 10)
    rowgap!(fig.layout, 1, 4)

    mkpath(FIGURES)
    record(fig, path, bearings; framerate = 12) do β
        compute!(β)
        label[] = SWEEP_LABEL[](β)
    end
    return path
end

"""
    figure_classes(config)

The same release under all six Pasquill classes, on one colour scale.
"""
function figure_classes(config)
    site = config.site
    panels = [
        ("$(letter(c))", (e, nn, β) -> dilution_instantaneous(e, nn, 0.0, site, c, β))
        for c in PASQUILL_CLASSES
    ]
    return panel_sweep(
        joinpath(FIGURES, "stability_classes.gif"),
        panels;
        columns = 3,
        width = 1080,
        height = 770,
        caption = "Pasquill class A (very unstable) to F (very stable), one colour scale",
    )
end

"""
    figure_buildings(config)

The same release with and without a building beside the stack, in the worst
case: one tall enough to swallow the plume entirely.

The wake criterion is a threshold, not a gradient. A stack clearing two and a
half building heights escapes untouched; one leaving *below* the building top is
entrained into the aerodynamic cavity and released at ground level. A 60 m
building next to this 50.3 m stack is on the far side of that threshold, so the
effective release height goes from 50.3 m to **zero** and the comparison shows
the whole effect rather than a fraction of it.
"""
function figure_buildings(config)
    source, air = config.site.source, config.site.atmosphere
    bare = Site(; source, atmosphere = air)
    east, north, h = 25.0, 0.0, 60.0
    building = Building(; east, north, height = h, frontal_area = 3600.0)
    waked = Site(; source, atmosphere = air, buildings = BuildingEnvelope([building]))
    @printf(
        "  building wake: release height %.1f m bare, %.1f m waked\n",
        release_height(bare),
        release_height(waked)
    )
    return comparison_sweep(
        joinpath(FIGURES, "building_wake.gif"),
        (
            "No building",
            (e, nn, β) -> dilution_instantaneous(e, nn, 0.0, bare, PASQUILL_D, β),
        ),
        (
            "With a 60 m building",
            (e, nn, β) -> dilution_instantaneous(e, nn, 0.0, waked, PASQUILL_D, β),
        );
        half_width = 12_000.0,
        caption = "Class D, 60 m building 25 m from a 50.3 m stack. The plume is entrained into the cavity and released at ground level: 50.3 m becomes 0 m",
    )
end

"""
    figure_heights(config)

The same release from three stack heights.
"""
function figure_heights(config)
    s, air = config.site.source, config.site.atmosphere
    sites = map((30.0, 50.3, 120.0)) do h
        source = StackSource(;
            height = h,
            diameter = s.diameter,
            exit_velocity = s.exit_velocity,
            exit_density = s.exit_density,
            exit_temperature = s.exit_temperature,
        )
        Site(; source, atmosphere = air)
    end
    panels = [
        (
            "$(Int(round(h))) m",
            (e, nn, β) -> dilution_instantaneous(e, nn, 0.0, st, PASQUILL_D, β),
        ) for (h, st) in zip((30.0, 50.3, 120.0), sites)
    ]
    return panel_sweep(
        joinpath(FIGURES, "release_height.gif"),
        panels;
        columns = 3,
        width = 1080,
        height = 450,
        caption = "Stack height, class D. A higher release moves the ground-level maximum downwind",
    )
end

"""
    figure_depletion(config)

The same release with and without depletion by decay, deposition and washout.

Depletion is a multiplicative factor below one that grows with distance, so on a
three-decade colour ramp the two fields look the same. The ratio panel is what
shows it.
"""
function figure_depletion(config)
    site, nuclide = config.site, config.nuclide
    washout = 3600.0
    depleted(e, nn, β) = dilution_instantaneous(
        e,
        nn,
        0.0,
        site,
        PASQUILL_D,
        β;
        nuclide,
        washout_duration = washout,
        precipitation = config.precipitation,
        rate = config.precipitation_rate,
    )
    return comparison_sweep(
        joinpath(FIGURES, "depletion.gif"),
        (
            "Undepleted",
            (e, nn, β) -> dilution_instantaneous(e, nn, 0.0, site, PASQUILL_D, β),
        ),
        ("Depleted", depleted);
        half_width = 20_000.0,
        ratio_span = 3.0,
        caption = "Class D, $(nuclide.name), one hour of washout at the configured rate",
    )
end

function main()
    config = reference_case()
    set_sweep_label!(config.rose)
    p1, sector = figure_field(config)
    println("wrote ", p1, "   (most exposed sector ", sector, ")")
    println("wrote ", figure_regimes(config))
    println("wrote ", figure_animation(config))
    println("wrote ", figure_classes(config))
    println("wrote ", figure_buildings(config))
    println("wrote ", figure_heights(config))
    println("wrote ", figure_depletion(config))
end

main()
