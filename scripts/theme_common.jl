# Figure style shared by every script in this directory, so the README figures
# and the validation figures are drawn to one standard.

using CairoMakie, MathTeXEngine

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
