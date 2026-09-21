# Figure style shared by every script in this directory, so that the README
# figures and the validation figures are drawn to one standard: fonts, sizing,
# axis treatment and the small set of formatting helpers below.

using CairoMakie, MathTeXEngine, LaTeXStrings

# One layout standard for every figure. A single panel is drawn on
# `Figure(size = (900, 600))`; each further stacked main panel adds about 350 to
# the height and each auxiliary strip (ratio, residual) about 180. Scripts scale
# up from these values where a figure needs it and never below them: no
# `fontsize`, `linewidth` or `markersize` override may undercut the theme,
# except `ANNOTATION_SIZE` for in-axis notes and the `cloud` and `dense` marker
# sizes below.
set_theme!(
    Theme(
    fonts = (;
        regular = texfont(:text),
        bold = texfont(:bold),
        italic = texfont(:italic),
    ),
    fontsize = 26,
    figure_padding = 10,
    linewidth = 3,
    markersize = 14,
    rowgap = 10,
    colgap = 12,
    Axis = (
        spinewidth = 1.5,
        xticklabelsize = 22,
        yticklabelsize = 22,
        xlabelpadding = 8,
        ylabelpadding = 8,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xgridcolor = (:grey, 0.12),
        ygridcolor = (:grey, 0.12),
        xminorticksvisible = false,
        yminorticksvisible = false,
        xtickalign = 1,
        ytickalign = 1,
    ),
    Scatter = (strokewidth = 1.5,),
    Legend = (
        framevisible = false,
        orientation = :horizontal,
        titlefont = :bold,
        labelsize = 22,
        padding = (0, 0, 0, 0),
    ),
    Colorbar = (ticklabelsize = 22, spinewidth = 1.5),
),
)

"Size of in-axis annotations, 0.8 of the base size. Nothing in a figure is set smaller."
const ANNOTATION_SIZE = 21

"Width of reference and guide lines, which are also dashed; data lines take the theme's 3."
const GUIDE_WIDTH = 1.5

"Okabe–Ito, colourblind-safe. One consistent colour per quantity across the figures."
#! format: off
const PALETTE = (
    blue   = "#0072B2",
    orange = "#E69F00",
    green  = "#009E73",
    red    = "#D55E00",
    purple = "#CC79A7",
    sky    = "#56B4E9",
    yellow = "#F0E442",
    black  = "#000000",
)
#! format: on

"""
Marker sizes, one scale for every figure, chosen against the shared
`fontsize = 26` so that symbols stay legible when a figure is reduced to
journal column width.

  * `cloud` — scatter of thousands of points, where larger markers would only
    deepen the overplotting
  * `dense` — series of a hundred points or more, where full-size markers would
    collide
  * `data` — the default for a measured or computed series
  * `emphasis` — single highlighted points: a fixed point, a fitted optimum
  * `key` — markers built by hand for a legend, which must stay readable even
    where the series itself is drawn small
"""
const MARKERSIZE = (cloud = 6, dense = 10, data = 14, emphasis = 20, key = 16)

_decimal_label(e::Integer) = e >= 0 ? latexstring(string(10^e)) :
                             latexstring("0." * "0"^(-e - 1) * "1")

# The unit decade and the first are written as the numbers they are: `10^0` and
# `10^1` are never the readable form, on an axis of powers no more than on one
# of decimals.
_power_label(e::Integer) = e == 0 ? L"1" : e == 1 ? L"10" : latexstring("10^{$e}")

"""
    logticks(e_min, e_max; step = 1, style = :auto)

Explicit decade ticks for a base-10 logarithmic axis running from `10^e_min` to
`10^e_max`, returned as the `(values, labels)` tuple an `Axis` takes.

Spans of four decades or fewer are labelled as plain decimals (`0.01, 0.1, 1,
10`) provided the decades themselves stay within `10^-4` to `10^4`, beyond which
a decimal is a run of zeros and the power of ten is the readable form; wider
spans use powers of ten throughout. Pass `style = :decimal` or `style = :power`
to override the choice.

Every logarithmic axis in these scripts sets its ticks through this function.
Left to itself Makie labels the unit decade as `10^0` on some axes and as `1` on
others, and on ranges narrower than a decade it falls back to fractional
exponents such as `10^{-0.4}`, which is not a form any of these quantities are
read in.

In the power style the unit decade is still written `1` and the first `10`, as
[`sci`](@ref) writes them: `10^{-2}, 1, 10^{2}`.

# Example

```julia
Axis(fig[1, 1]; xscale = log10, xticks = logticks(-4, 0))
```
"""
function logticks(e_min::Integer, e_max::Integer; step::Integer = 1, style::Symbol = :auto)
    e_min <= e_max || throw(ArgumentError("e_min = $e_min exceeds e_max = $e_max"))
    step > 0 || throw(ArgumentError("step must be positive, got $step"))
    style in (:auto, :decimal, :power) ||
        throw(ArgumentError("style must be :auto, :decimal or :power, got :$style"))
    exponents = collect(e_min:step:e_max)
    readable = e_max - e_min <= 4 && e_min >= -4 && e_max <= 4
    decimal = style == :decimal || (style == :auto && readable)
    labels = [decimal ? _decimal_label(e) : _power_label(e) for e in exponents]
    return ([10.0^e for e in exponents], labels)
end

"""
    it(s)

`s` set in the italic of the shared Computer Modern family, for use inside a
`rich(...)` run: `rich("Lorenz: ", it("σ"), " = 10")`.

`rich` rather than a `LaTeXString` wherever a label carries ordinary prose,
because MathTeXEngine parses everything in an `L"…"` — including the text
outside `\$…\$` — as mathematics, so a colon picks up relation spacing and
renders as `Lorenz : σ = 10`. Reserve `L"…"` for labels that are entirely
mathematical.
"""
it(s) = rich(s, font = :italic)

_trim(v::Real) = (s = string(v); endswith(s, ".0") ? s[1:(end - 2)] : s)

"""
    sci(x; digits = 2)

`x` written as a LaTeX power of ten — `sci(1.8e-7)` gives `1.8 \\times 10^{-7}`.

Exponents of −1, 0 and 1 are folded into the mantissa (`0.18`, `1.8`, `18`) and a
mantissa of 1 is dropped, so the result never reads `1 \\times 10^{n}` or carries
a `\\times 10^{0}`. Figures and papers never show computer notation such as
`1.8e-07`, which is what `@sprintf("%e", …)` produces.
"""
function sci end

"""
    _split_sci(x, digits) -> (mantissa::String, exponent::Int) or (folded::String, nothing)

The mantissa and exponent `sci` and `rsci` share. Exponents of −1, 0 and 1 come
back already folded into a plain decimal, with `nothing` for the exponent.
"""
function _split_sci(x::Real, digits::Integer)
    isfinite(x) || throw(ArgumentError("expected a finite number, got $x"))
    e = floor(Int, log10(abs(x)))
    m = round(x / 10.0^e; digits = digits)
    if abs(m) >= 10                      # rounding can carry the mantissa to 10
        m /= 10
        e += 1
    end
    -1 <= e <= 1 && return (_trim(round(x; sigdigits = digits + 1)), nothing)
    return (_trim(m), e)
end

function sci(x::Real; digits::Integer = 2)
    iszero(x) && return L"0"
    mantissa, e = _split_sci(x, digits)
    e === nothing && return latexstring(mantissa)
    mantissa == "1" && return latexstring("10^{$e}")
    mantissa == "-1" && return latexstring("-10^{$e}")
    return latexstring("$mantissa \\times 10^{$e}")
end

"""
    rsci(x; digits = 2)

`sci` as a `rich` run, for the mixed prose-and-mathematics labels a
`LaTeXString` cannot be embedded in. Minus signs are the typographic `−`.
"""
function rsci(x::Real; digits::Integer = 2)
    iszero(x) && return rich("0")
    mantissa, e = _split_sci(x, digits)
    minus(s) = replace(s, "-" => "−")
    e === nothing && return rich(minus(mantissa))
    mantissa == "1" && return rich("10", superscript(minus(string(e))))
    mantissa == "-1" && return rich("−10", superscript(minus(string(e))))
    return rich(minus(mantissa), " × 10", superscript(minus(string(e))))
end

"""
    savefigure(fig, dir, name)

Write `fig` to `dir/name.png` well above 300 dpi, creating `dir` if needed.
Callers pass `joinpath(@__DIR__, "figures")` so that nothing depends on the
working directory. Returns the path written.
"""
function savefigure(fig, dir::AbstractString, name::AbstractString)
    mkpath(dir)
    path = joinpath(dir, name * ".png")
    save(path, fig; px_per_unit = 4)
    return path
end
