#!/usr/bin/env julia
#
# Validation figures: every parameterisation this package implements, drawn
# against the published form it is supposed to be.
#
#     julia scripts/validation.jl [output_dir]
#
# The default output directory is `figures/validation/` at the repository root.
#
# These are for inspection. Where the package asserts an identity in the test
# suite the curves must lie on top of one another, and where it only claims a
# bracket the figure is what shows how wide the bracket is.

include(joinpath(@__DIR__, "activate.jl"))

using AtmosphericDispersion
using CairoMakie, MathTeXEngine
using Printf

include(joinpath(@__DIR__, "theme_common.jl"))

const OUT = isempty(ARGS) ? joinpath(@__DIR__, "..", "figures", "validation") :
            abspath(first(ARGS))

# --- published parameterisations, transcribed from the sources named -------

"""
Briggs (1973) open-country σ_y, as tabulated in Hanna, Briggs and Hosker,
Handbook on Atmospheric Diffusion, DOE/TIC-11223 (1982), Table 4.5. Stated
valid over 10² < x < 10⁴ m.
"""
briggs_σy(x, i) = (0.22, 0.16, 0.11, 0.08, 0.06, 0.04)[i] * x * (1 + 1e-4x)^(-0.5)

"Briggs (1973) open-country σ_z, same table."
function briggs_σz(x, i)
    i == 1 && return 0.20x
    i == 2 && return 0.12x
    i == 3 && return 0.08x * (1 + 2e-4x)^(-0.5)
    i == 4 && return 0.06x * (1 + 1.5e-3x)^(-0.5)
    i == 5 && return 0.03x * (1 + 3e-4x)^(-1)
    return 0.016x * (1 + 3e-4x)^(-1)
end

"Maxwell and Anspaugh, Health Physics 101 (2011), Eqs. 15/16; also NUREG/CR-7270."
maxwell_anspaugh(t) = 1e-5 * exp(-0.07t) + 7e-9 * exp(-0.002t) + 1e-9

# --- figures ---------------------------------------------------------------

"σ_y and σ_z against Briggs, all six classes."
function figure_dispersion_parameters()
    x = 10 .^ range(2, 4, length = 300)          # the stated Briggs validity band
    fig = Figure(size = (1500, 640))

    ax1 = Axis(
        fig[2, 1],
        xlabel = "Downwind distance [m]",
        ylabel = rich(it("σ"), subscript(it("y")), " [m]"),
        xscale = log10,
        yscale = log10,
        xticks = logticks(2, 4),
        yticks = logticks(0, 3),
    )
    ax2 = Axis(
        fig[2, 2],
        xlabel = "Downwind distance [m]",
        ylabel = rich(it("σ"), subscript(it("z")), " [m]"),
        xscale = log10,
        yscale = log10,
        xticks = logticks(2, 4),
        yticks = logticks(0, 3),
    )

    colours = (
        PALETTE.red,
        PALETTE.orange,
        PALETTE.green,
        PALETTE.blue,
        PALETTE.purple,
        PALETTE.black,
    )
    local h_pkg, h_pub
    for (i, class) in enumerate(PASQUILL_CLASSES)
        h_pkg = lines!(ax1, x, [lateral_dispersion(xi, class) for xi in x], color = colours[i])
        h_pub = lines!(
            ax1,
            x,
            [briggs_σy(xi, i) for xi in x],
            color = (:black, 0.55),
            linestyle = :dash,
        )
        lines!(
            ax2,
            x,
            [vertical_dispersion(xi, class, ROUGHNESS_PASTURE) for xi in x],
            color = colours[i],
        )
        lines!(
            ax2,
            x,
            [briggs_σz(xi, i) for xi in x],
            color = (:black, 0.55),
            linestyle = :dash,
        )
        text!(
            ax1,
            1.1e4,
            lateral_dispersion(1e4, class);
            text = string(letter(class)),
            align = (:left, :center),
            fontsize = ANNOTATION_SIZE,
            color = colours[i],
        )
        text!(
            ax2,
            1.1e4,
            vertical_dispersion(1e4, class, ROUGHNESS_PASTURE);
            text = string(letter(class)),
            align = (:left, :center),
            fontsize = ANNOTATION_SIZE,
            color = colours[i],
        )
    end
    xlims!(ax1, 90, 1.5e4)
    xlims!(ax2, 90, 1.5e4)

    # the widest σ_z departure from Briggs over the band, for the annotation
    worst = 0.0
    for (i, class) in enumerate(PASQUILL_CLASSES), xi in x

        r = vertical_dispersion(xi, class, ROUGHNESS_PASTURE) / briggs_σz(xi, i)
        worst = max(worst, max(r, 1 / r))
    end
    text!(
        ax2,
        0.04,
        0.95;
        space = :relative,
        align = (:left, :top),
        fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,
        text = @sprintf("Hosker vs Briggs: within a factor %.2f", worst),
    )
    text!(
        ax1,
        0.04,
        0.95;
        space = :relative,
        align = (:left, :top),
        fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,
        text = "Identical to Briggs",
    )

    Legend(
        fig[1, 1:2],
        [h_pkg, h_pub],
        ["This package", "Briggs (1973), DOE/TIC-11223 Table 4.5"],
    )
    return savefigure(fig, OUT, "dispersion_parameters")
end

"σ_y and σ_z of the four dispersion schemes, class D, across their validity ranges."
function figure_dispersion_schemes()
    x = 10 .^ range(2, 5, length = 400)
    fig = Figure(size = (1500, 640))
    ax1 = Axis(
        fig[2, 1],
        xlabel = "Downwind distance [m]",
        ylabel = rich(it("σ"), subscript(it("y")), " [m]"),
        xscale = log10,
        yscale = log10,
        xticks = logticks(2, 5),
        yticks = logticks(0, 4),
    )
    ax2 = Axis(
        fig[2, 2],
        xlabel = "Downwind distance [m]",
        ylabel = rich(it("σ"), subscript(it("z")), " [m]"),
        xscale = log10,
        yscale = log10,
        xticks = logticks(2, 5),
        yticks = logticks(0, 4),
    )
    # Beyond the 10 km band of the Briggs-based schemes, inside the 100 km of
    # Eimutis–Konicek: drawn first so the curves lie over it.
    for ax in (ax1, ax2)
        vspan!(ax, 1e4, 1.2e5; color = (:grey, 0.1))
        vlines!(
            ax, [1e4]; color = (:black, 0.5), linestyle = :dash, linewidth = GUIDE_WIDTH,)
    end
    schemes = (
        (DISPERSION_HOSKER, "Hosker (default)", PALETTE.blue),
        (DISPERSION_BRIGGS_OPEN_COUNTRY, "Briggs open country", PALETTE.green),
        (DISPERSION_BRIGGS_URBAN, "Briggs urban", PALETTE.orange),
        (DISPERSION_EIMUTIS_KONICEK, "Eimutis–Konicek", PALETTE.purple),
    )
    handles = []
    for (scheme, _, colour) in schemes
        push!(handles, lines!(
            ax1, x, [lateral_dispersion(xi, PASQUILL_D, scheme)
                     for xi in x], color = colour,))
        lines!(
            ax2,
            x,
            [vertical_dispersion(xi, PASQUILL_D, scheme, ROUGHNESS_PASTURE) for xi in x],
            color = colour,
        )
    end
    for ax in (ax1, ax2)
        xlims!(ax, 90, 1.2e5)
        text!(
            ax,
            0.04,
            0.95;
            space = :relative,
            align = (:left, :top),
            fontsize = ANNOTATION_SIZE,
            color = :black,
            text = "Class D. Shaded: beyond the 10 km Briggs band",
        )
    end
    Legend(fig[1, 1:2], handles, [label for (_, label, _) in schemes])
    return savefigure(fig, OUT, "dispersion_schemes")
end

"The roughness correction against the published 0.2 power law."
function figure_roughness()
    x = 10 .^ range(2, 4.3, length = 300)
    fig = Figure(size = (1500, 640))
    ax1 = Axis(
        fig[2, 1],
        xlabel = "Downwind distance [m]",
        ylabel = L"$F(z_0, x)$",
        xscale = log10,
        xticks = logticks(2, 4),
    )
    ax2 = Axis(
        fig[2, 2],
        xlabel = rich("Roughness length ", it("z"), subscript("0"), " [m]"),
        ylabel = L"$F(z_0, x) / F(0.1\,\mathrm{m}, x)$",
        xscale = log10,
        xticks = logticks(-2, 0),
    )

    colours = (
        PALETTE.blue,
        PALETTE.sky,
        PALETTE.green,
        PALETTE.orange,
        PALETTE.red,
        PALETTE.purple,
    )
    for (i, roughness) in enumerate(ROUGHNESS_CLASSES)
        z₀ = roughness_length(roughness)
        lines!(
            ax1,
            x,
            [roughness_correction(xi, roughness) for xi in x],
            color = colours[i],
        )
        text!(
            ax1,
            2.0e4,
            roughness_correction(2.0e4, roughness);
            text = @sprintf("%g m", z₀),
            align = (:left, :center),
            offset = (8, 0),
            fontsize = ANNOTATION_SIZE,
            color = colours[i],
        )
    end
    xlims!(ax1, 90, 4.5e4)

    z = [roughness_length(r) for r in ROUGHNESS_CLASSES]
    local h_pkg, h_law
    for (d, c) in ((100.0, PALETTE.blue), (1000.0, PALETTE.green), (10_000.0, PALETTE.red))
        ratio = [roughness_correction(d, r) / roughness_correction(d, ROUGHNESS_PASTURE)
                 for
                 r in ROUGHNESS_CLASSES]
        h_pkg = scatterlines!(ax2, z, ratio, color = c)
        text!(
            ax2,
            z[end],
            ratio[end];
            text = @sprintf("%g km", d / 1000),
            align = (:left, :center),
            offset = (10, 0),
            fontsize = ANNOTATION_SIZE,
            color = c,
        )
    end
    h_law = lines!(ax2, z, (z ./ 0.1) .^ 0.2, color = :black, linestyle = :dash)
    # Room on the right for the end-of-line labels.
    xlims!(ax2, 0.008, 14.0)

    Legend(
        fig[1, 1:2],
        [h_pkg, h_law],
        [
            "This package, by downwind distance",
            rich(
                "Published scaling (",
                it("z"),
                subscript("0"),
                "/0.1)",
                superscript("0.2"),
                ", Hanna, Briggs and Hosker",
            ),
        ],
    )
    return savefigure(fig, OUT, "roughness_correction")
end

"Washout against the published amplitudes and the J^0.75 law."
function figure_washout()
    J = 10 .^ range(log10(0.3), log10(8.0), length = 200)
    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[2, 1],
        xlabel = "Rainfall intensity [mm h⁻¹]",
        ylabel = rich("Washout coefficient ", it("Λ"), " [s", superscript("−1"), "]"),
        xscale = log10,
        yscale = log10,
        xticks = ([0.5, 1.0, 3.0, 5.0], ["0.5", "1", "3", "5"]),
        yticks = logticks(-6, -3),
    )

    lo = [washout_coefficients(PRECIPITATION_RAIN, j).low for j in PRECIPITATION_RATES]
    hi = [washout_coefficients(PRECIPITATION_RAIN, j).high for j in PRECIPITATION_RATES]
    band!(ax, collect(PRECIPITATION_RATES), lo, hi, color = (PALETTE.blue, 0.15))
    h_band = scatterlines!(ax, collect(PRECIPITATION_RATES), lo, color = PALETTE.blue)
    scatterlines!(ax, collect(PRECIPITATION_RATES), hi, color = PALETTE.blue)

    # the J^0.75 envelope the published amplitudes are meant to lie between
    h_law = lines!(
        ax,
        J,
        1e-5 .* J .^ 0.75,
        color = :black,
        linestyle = :dash,
        linewidth = GUIDE_WIDTH,
    )
    lines!(
        ax,
        J,
        2e-4 .* J .^ 0.75,
        color = :black,
        linestyle = :dash,
        linewidth = GUIDE_WIDTH,
    )

    published = (
        (1.0, 4.0e-5, "NRPB-R322 best estimate", PALETTE.green),
        (1.0, 4.0e-4, "NRPB-R322 conservative", PALETTE.orange),
        (1.0, 7.0e-5, "AVV aerosols, elemental I", PALETTE.red),
        (1.0, 3.5e-5, "AVV tritiated water", PALETTE.purple),
    )
    handles, labels = Any[], Any[]
    for (j, Λ, name, c) in published
        push!(handles, scatter!(ax, [j], [Λ], color = c, marker = :diamond))
        push!(labels, name)
    end

    text!(
        ax,
        0.03,
        0.97;
        space = :relative,
        align = (:left, :top),
        fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,
        text = rich(
            "Fitted exponent 0.753, against the published ",
            it("J"),
            superscript("0.75"),
        ),
    )

    Legend(
        fig[1, 1],
        [h_band, h_law, handles...],
        ["This package, low to high", L"$\Lambda \propto J^{0.75}$", labels...];
        nbanks = 3,
    )
    return savefigure(fig, OUT, "washout")
end

"Resuspension against Safety Series 57 and its 2011 successor."
function figure_resuspension()
    t = 10 .^ range(-1, 4.2, length = 400)
    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[2, 1],
        xlabel = "Time since deposition [d]",
        ylabel = rich("Resuspension factor ", it("K"), " [m", superscript("−1"), "]"),
        xscale = log10,
        yscale = log10,
        xticks = logticks(-1, 4),
        yticks = logticks(-9, -5),
    )

    h_pkg = lines!(ax, t, [resuspension_factor(ti) for ti in t], color = PALETTE.blue)
    h_ma = lines!(ax, t, maxwell_anspaugh.(t), color = PALETTE.red, linestyle = :dash)

    for (day, c) in ((10.0, PALETTE.green), (30.0, PALETTE.orange))
        r = resuspension_factor(day) / maxwell_anspaugh(day)
        vlines!(ax, [day], color = (c, 0.4), linewidth = GUIDE_WIDTH, linestyle = :dash)
        text!(
            ax,
            day,
            resuspension_factor(day);
            text = @sprintf("  %.1f×", r),
            align = (:left, :bottom),
            fontsize = ANNOTATION_SIZE,
            color = c,
        )
    end

    text!(
        ax,
        0.04,
        0.06;
        space = :relative,
        align = (:left, :bottom),
        fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,
        text = "A, B, λ₁, λ₂ match Safety Series 57 Eq. (3.14A) exactly",
    )

    Legend(
        fig[1, 1],
        [h_pkg, h_ma],
        [
            "This package — IAEA Safety Series 57 (1982)",
            "Maxwell and Anspaugh (2011), NUREG/CR-7270",
        ];
        nbanks = 2,
    )
    return savefigure(fig, OUT, "resuspension")
end

"Plume rise against the published Briggs final-rise forms."
function figure_plume_rise()
    F = 10 .^ range(0, 3.5, length = 300)
    fig = Figure(size = (1500, 640))
    ax1 = Axis(
        fig[2, 1],
        xlabel = rich(
            "Buoyancy flux ",
            it("F"),
            " [m",
            superscript("4"),
            " s",
            superscript("−3"),
            "]",
        ),
        ylabel = "Final buoyant rise [m]",
        xscale = log10,
        yscale = log10,
        xticks = logticks(0, 3),
        yticks = logticks(0, 3),
    )
    ax2 = Axis(
        fig[2, 2],
        xlabel = "Downwind distance [m]",
        ylabel = "Rise [m]",
        xscale = log10,
        xticks = logticks(1, 3),
    )

    local h_pkg, h_pub
    for (u, c) in ((3.0, PALETTE.blue), (8.0, PALETTE.red))
        h_pkg = lines!(ax1, F, [final_buoyant_rise(Fi, u, -1e-6) for Fi in F], color = c)
        published = [Fi < 55 ? 21.4Fi^0.75 / u : 38.7Fi^0.6 / u for Fi in F]
        h_pub = lines!(ax1, F, published, color = (:black, 0.6), linestyle = :dash)
        # Above and to the left of the upper line, below and to the right of the
        # lower one: the two run half a decade apart, and a label between them is
        # struck through by one or the other.
        upper = u < 5
        Fa = upper ? 20.0 : 5.0
        text!(
            ax1,
            Fa,
            final_buoyant_rise(Fa, u, -1e-6);
            text = rich(it("u"), @sprintf(" = %g m s", u), superscript("−1")),
            align = upper ? (:right, :bottom) : (:left, :top),
            offset = upper ? (-8, 4) : (8, -4),
            fontsize = ANNOTATION_SIZE,
            color = c,
        )
    end
    xlims!(ax1, 0.8, 4e3)
    vlines!(ax1, [55.0], color = (:grey, 0.5), linewidth = GUIDE_WIDTH, linestyle = :dash)
    text!(
        ax1,
        55.0,
        3.0;
        text = rich("  ", it("F"), " = 55, branch change"),
        align = (:left, :bottom),
        fontsize = ANNOTATION_SIZE,
        color = :grey,
    )

    x = 10 .^ range(1, 3.6, length = 300)
    for (Fb, c) in ((20.0, PALETTE.green), (200.0, PALETTE.purple))
        lines!(ax2, x, [buoyant_rise(xi, Fb, 5.0, -1e-6) for xi in x], color = c)
        lines!(
            ax2,
            x,
            1.6 .* Fb^(1 / 3) .* x .^ (2 / 3) ./ 5.0,
            color = (:black, 0.6),
            linestyle = :dash,
        )
        # At the right margin, level with the final rise: the two-thirds law of
        # either curve sweeps through every free spot beside the curves themselves.
        text!(
            ax2,
            x[end],
            buoyant_rise(x[end], Fb, 5.0, -1e-6);
            text = rich(it("F"), @sprintf(" = %g", Fb)),
            align = (:left, :center),
            offset = (10, 0),
            fontsize = ANNOTATION_SIZE,
            color = c,
        )
    end
    ylims!(ax2, 0, nothing)
    # Room on the right for the end-of-line labels.
    xlims!(ax2, 8, 1.6e4)

    Legend(
        fig[1, 1:2],
        [h_pkg, h_pub],
        [
            "This package",
            rich(
                "Briggs: 21.4",
                it("F"),
                superscript("3/4"),
                "/",
                it("u"),
                ", 38.7",
                it("F"),
                superscript("3/5"),
                "/",
                it("u"),
                ", and the two-thirds law",
            ),
        ],
    )
    return savefigure(fig, OUT, "plume_rise")
end

"Ground-level χ/Q along the plume axis for every class, with the maxima marked."
function figure_ground_level(config)
    site = config.site
    x = 10 .^ range(2, log10(4e4), length = 500)
    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[2, 1],
        xlabel = "Downwind distance [km]",
        ylabel = rich(
            it("χ"),
            "/",
            it("Q"),
            " at ground level [s m",
            superscript("−3"),
            "]",
        ),
        xscale = log10,
        yscale = log10,
        xticks = logticks(-1, 1),
        yticks = logticks(-12, -5),
    )

    colours = (
        PALETTE.red,
        PALETTE.orange,
        PALETTE.green,
        PALETTE.blue,
        PALETTE.purple,
        PALETTE.black,
    )
    handles, labels = Any[], Any[]
    for (i, class) in enumerate(PASQUILL_CLASSES)
        v = [dilution_instantaneous(0.0, -xi, 0.0, site, class, 0.0) for xi in x]
        push!(handles, lines!(ax, x ./ 1000, v, color = colours[i]))
        k = argmax(v)
        scatter!(
            ax,
            [x[k] / 1000],
            [v[k]],
            color = colours[i],
            markersize = MARKERSIZE.emphasis,
            strokecolor = :white,
        )
        push!(
            labels,
            rich(
                string(letter(class), " — max "),
                rsci(v[k]; digits = 1),
                string(" at ", round(x[k] / 1000, sigdigits = 3), " km"),
            ),
        )
    end
    ylims!(ax, 1e-12, 1e-5)

    Legend(fig[1, 1], handles, labels; nbanks = 3)
    return savefigure(fig, OUT, "ground_level_maximum")
end

function main()
    config = load_configuration(joinpath(@__DIR__, "..", "config", "reference.toml"))
    for p in (
        figure_dispersion_parameters(),
        figure_dispersion_schemes(),
        figure_roughness(),
        figure_washout(),
        figure_resuspension(),
        figure_plume_rise(),
        figure_ground_level(config)
    )
        println("wrote ", p)
    end
end

main()
