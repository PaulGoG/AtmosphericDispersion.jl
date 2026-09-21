#=
Configuration.

A run is described by a TOML file, not by editing source. The loader enforces
exactly what the file's comments promise — presence, type, enumerated choice,
numerical bound, and that no key is present which the loader does not read —
and every failure names the key it failed on, by its full dotted path, so that
a rejected file says what to change rather than where the parser happened to
stop. An unread key is an error because it is otherwise indistinguishable from
a misspelt one, which would fall back to its default without a word.

A pipeline must be unable to start from a configuration it cannot honour, so
the checks run at load and construct the solver types immediately: anything the
constructors themselves reject (a negative stack height, frequencies that do
not sum to one) fails here too.
=#

"""
    ConfigurationError(path, message)

A configuration file that could not be honoured. `path` is the dotted key the
failure is attached to; `message` says what was wrong with it.
"""
struct ConfigurationError <: Exception
    path::String
    message::String
end

Base.showerror(io::IO, e::ConfigurationError) = print(
    io, "ConfigurationError at `", e.path, "`: ", e.message,)

_fail(path, message) = throw(ConfigurationError(path, message))

function _table(parent::AbstractDict, key::AbstractString, path::AbstractString)
    haskey(parent, key) || _fail(_join(path, key), "required table is missing")
    value = parent[key]
    value isa AbstractDict ||
        _fail(_join(path, key), "expected a table, got $(typeof(value))")
    return value
end

_join(path, key) = isempty(path) ? String(key) : string(path, ".", key)

function _value(
        parent::AbstractDict,
        key::AbstractString,
        ::Type{T},
        path::AbstractString;
        default = nothing,
) where {T}
    if !haskey(parent, key)
        default === nothing && _fail(_join(path, key), "required key is missing")
        return default::T
    end
    value = parent[key]
    value isa T && return value
    # TOML gives integers for unadorned numerals; accept them where a float is wanted.
    T === Float64 && value isa Integer && return Float64(value)
    return _fail(_join(path, key), "expected $T, got $(typeof(value))")
end

function _positive(value::Real, path::AbstractString)
    value > 0 || _fail(path, "must be positive, got $value")
    return value
end

function _nonnegative(value::Real, path::AbstractString)
    value ≥ 0 || _fail(path, "cannot be negative, got $value")
    return value
end

function _choice(value::AbstractString, options::AbstractDict, path::AbstractString)
    haskey(options, value) && return options[value]
    choices = join(sort(collect(keys(options))), ", ")
    return _fail(path, "must be one of $choices, got $(repr(value))")
end

function _reject_unknown(t::AbstractDict, allowed, path::AbstractString)
    for key in sort!(collect(keys(t)))
        key in allowed && continue
        _fail(_join(path, key), "unknown key; expected one of $(join(allowed, ", "))")
    end
    return nothing
end

const _ROOT_KEYS = (
    "source",
    "atmosphere",
    "buildings",
    "model",
    "mixing_layer",
    "nuclide",
    "wind_rose",
    "release",
    "precipitation",
    "grid",
)
const _SOURCE_KEYS = (
    "height", "diameter", "exit_velocity", "exit_density", "exit_temperature",)
const _ATMOSPHERE_KEYS = (
    "reference_speed",
    "temperature",
    "density",
    "lapse_rate",
    "specific_heat",
    "surface",
    "roughness",
)
const _MODEL_KEYS = ("plume_rise", "resuspension", "washout")
const _MIXING_KEYS = ("scheme", "above_lid", "uniform_depth", "depths")
const _BUILDINGS_KEYS = ("wake_coefficient", "building")
const _BUILDING_KEYS = ("east", "north", "height", "frontal_area")
const _NUCLIDE_KEYS = (
    "name",
    "decay_constant",
    "deposition_velocity_low",
    "deposition_velocity_high",
    "washout_species",
)
const _ROSE_KEYS = ("sectors", "convention", "frequencies", "stability")
const _RELEASE_KEYS = ("activity", "duration")
const _PRECIPITATION_KEYS = ("type", "rate", "washout_duration")
const _GRID_KEYS = ("extent", "spacing")

const _RISE_CHOICES = Dict("briggs" => BRIGGS_RISE, "xoqdoq" => XOQDOQ_RISE, "nsr23" =>
    NSR23_RISE)

const _MIXING_SCHEMES = ("tabulated", "uniform", "custom", "unbounded")

const _LID_CHOICES = Dict("rise_inhibited" => RISE_INHIBITED, "full_penetration" =>
    FULL_PENETRATION)

const _WASHOUT_CHOICES = Dict("normative" => WASHOUT_NORMATIVE, "hto" => WASHOUT_HTO)

const _SPECIES_CHOICES = Dict("tritium_iodine" => WASHOUT_TRITIUM_IODINE, "other" =>
    WASHOUT_OTHER_NUCLIDES)

const _RESUSPENSION_CHOICES = Dict(
    "iaea_ss57" => RESUSPENSION_IAEA_SS57,
    "maxwell_anspaugh" => RESUSPENSION_MAXWELL_ANSPAUGH,
)

const _SURFACE_CHOICES = Dict(
    "water" => SURFACE_WATER,
    "agricultural" => SURFACE_AGRICULTURAL,
    "forest_urban" => SURFACE_FOREST_URBAN,
)

const _ROUGHNESS_CHOICES = Dict(
    "grassland_water" => ROUGHNESS_GRASSLAND_WATER,
    "arable" => ROUGHNESS_ARABLE,
    "pasture" => ROUGHNESS_PASTURE,
    "rural" => ROUGHNESS_RURAL,
    "forest_urban" => ROUGHNESS_FOREST_URBAN,
    "metropolis" => ROUGHNESS_METROPOLIS,
)

const _PRECIPITATION_CHOICES = Dict("rain" => PRECIPITATION_RAIN, "snow" =>
    PRECIPITATION_SNOW)

const _CONVENTION_CHOICES = Dict("blowing_from" => BlowingFrom(), "blowing_toward" =>
    BlowingToward())

"""
    RunConfiguration

Everything a dispersion run needs, assembled and validated from a TOML file by
[`load_configuration`](@ref).

- `site` — the [`Site`](@ref), with its source, atmosphere and buildings
- `rose` — the [`WindRose`](@ref), already resolved to blowing-towards
- `nuclide` — the released [`Nuclide`](@ref)
- `resuspension` — the [`ResuspensionModel`](@ref) to apply
- `washout` — the [`WashoutEvent`](@ref), of zero duration for a dry plume
- `activity` — released activity, Bq
- `release_duration` — duration of the release, s
- `extent`, `spacing` — the receptor grid, m
"""
struct RunConfiguration
    site::Site
    rose::WindRose{Float64}
    nuclide::Nuclide
    resuspension::ResuspensionModel
    washout::WashoutEvent
    activity::Float64
    release_duration::Float64
    extent::Float64
    spacing::Float64
end

"""
    load_configuration(path)

Read and validate a run configuration from the TOML file at `path`.

Throws a [`ConfigurationError`](@ref) naming the offending key on the first
constraint that fails, and an `ArgumentError` from the solver constructors for
anything they reject in turn.
"""
function load_configuration(path::AbstractString)
    isfile(path) || throw(ArgumentError("no configuration file at $path"))
    return configuration_from(TOML.parsefile(String(path)))
end

"""
    configuration_from(table)

Validate an already-parsed TOML table into a [`RunConfiguration`](@ref).
"""
function configuration_from(root::AbstractDict)
    _reject_unknown(root, _ROOT_KEYS, "")
    source = _source_from(_table(root, "source", ""))
    atmosphere = _atmosphere_from(_table(root, "atmosphere", ""))
    envelope = get(root, "buildings", Dict{String,Any}())
    envelope isa AbstractDict ||
        _fail("buildings", "expected a table, got $(typeof(envelope))")
    buildings = _buildings_from(envelope)

    model = get(root, "model", Dict{String,Any}())
    model isa AbstractDict || _fail("model", "expected a table, got $(typeof(model))")
    _reject_unknown(model, _MODEL_KEYS, "model")
    rise = _choice(
        _value(model, "plume_rise", String, "model"; default = "briggs"),
        _RISE_CHOICES,
        "model.plume_rise",
    )
    resuspension = _choice(
        _value(model, "resuspension", String, "model"; default = "iaea_ss57"),
        _RESUSPENSION_CHOICES,
        "model.resuspension",
    )
    washout_model = _choice(
        _value(model, "washout", String, "model"; default = "normative"),
        _WASHOUT_CHOICES,
        "model.washout",
    )

    lid = get(root, "mixing_layer", Dict{String,Any}())
    lid isa AbstractDict || _fail("mixing_layer", "expected a table, got $(typeof(lid))")
    mixing = _mixing_from(lid)

    site = Site(; source, atmosphere, buildings, rise, mixing)

    rose = _rose_from(_table(root, "wind_rose", ""))
    nuclide = _nuclide_from(_table(root, "nuclide", ""))

    release = _table(root, "release", "")
    _reject_unknown(release, _RELEASE_KEYS, "release")
    activity = _nonnegative(_value(release, "activity", Float64, "release"), "release.activity")
    release_duration = _positive(_value(release, "duration", Float64, "release"), "release.duration")

    precip = get(root, "precipitation", Dict{String,Any}())
    precip isa AbstractDict ||
        _fail("precipitation", "expected a table, got $(typeof(precip))")
    _reject_unknown(precip, _PRECIPITATION_KEYS, "precipitation")
    precipitation = _choice(
        _value(precip, "type", String, "precipitation"; default = "rain"),
        _PRECIPITATION_CHOICES,
        "precipitation.type",
    )
    rate = _value(precip, "rate", Float64, "precipitation"; default = 1.0)
    if !(rate in PRECIPITATION_RATES)
        rates = join(PRECIPITATION_RATES, ", ")
        _fail(
            "precipitation.rate",
            "the washout table is defined at $rates mm/h, got $rate",
        )
    end
    washout = _nonnegative(
        _value(precip, "washout_duration", Float64, "precipitation"; default = 0.0),
        "precipitation.washout_duration",
    )

    grid = _table(root, "grid", "")
    _reject_unknown(grid, _GRID_KEYS, "grid")
    extent = _positive(_value(grid, "extent", Float64, "grid"), "grid.extent")
    spacing = _positive(_value(grid, "spacing", Float64, "grid"), "grid.spacing")
    spacing < extent ||
        _fail("grid.spacing", "must be smaller than grid.extent ($extent m), got $spacing")

    event = WashoutEvent(; duration = washout, precipitation, rate, model = washout_model)

    return RunConfiguration(
        site,
        rose,
        nuclide,
        resuspension,
        event,
        activity,
        release_duration,
        extent,
        spacing,
    )
end

function _mixing_from(t::AbstractDict)
    _reject_unknown(t, _MIXING_KEYS, "mixing_layer")
    scheme = _value(t, "scheme", String, "mixing_layer"; default = "tabulated")
    scheme in _MIXING_SCHEMES || _fail(
        "mixing_layer.scheme",
        "must be one of $(join(_MIXING_SCHEMES, ", ")), got $(repr(scheme))",
    )
    above_lid = _choice(
        _value(t, "above_lid", String, "mixing_layer"; default = "rise_inhibited"),
        _LID_CHOICES,
        "mixing_layer.above_lid",
    )
    # A key the chosen scheme does not read is refused, like any other unread key.
    for (key, owner) in (("uniform_depth", "uniform"), ("depths", "custom"))
        haskey(t, key) &&
            scheme != owner &&
            _fail("mixing_layer.$key", "is read only with scheme = $(repr(owner))")
    end
    scheme == "tabulated" && return MixingLayer(MIXING_TABULATED.depths; above_lid)
    scheme == "unbounded" && return MixingLayer(Inf; above_lid)
    if scheme == "uniform"
        depth = _positive(
            _value(t, "uniform_depth", Float64, "mixing_layer"),
            "mixing_layer.uniform_depth",
        )
        return MixingLayer(depth; above_lid)
    end
    depths = _floatvector(t, "depths", "mixing_layer")
    nclasses = length(PASQUILL_CLASSES)
    length(depths) == nclasses || _fail(
        "mixing_layer.depths",
        "must hold one depth per Pasquill class ($nclasses), got $(length(depths))",
    )
    for (i, d) in enumerate(depths)
        d > 0 ||
            _fail("mixing_layer.depths[$i]", "must be positive, `inf` for none, got $d")
    end
    return MixingLayer(depths; above_lid)
end

function _source_from(t::AbstractDict)
    _reject_unknown(t, _SOURCE_KEYS, "source")
    return StackSource(;
        height = _positive(_value(t, "height", Float64, "source"), "source.height"),
        diameter = _positive(_value(t, "diameter", Float64, "source"), "source.diameter"),
        exit_velocity = _nonnegative(
            _value(t, "exit_velocity", Float64, "source"),
            "source.exit_velocity",
        ),
        exit_density = _positive(
            _value(t, "exit_density", Float64, "source"),
            "source.exit_density",
        ),
        exit_temperature = _positive(
            _value(t, "exit_temperature", Float64, "source"),
            "source.exit_temperature",
        ),
    )
end

function _atmosphere_from(t::AbstractDict)
    _reject_unknown(t, _ATMOSPHERE_KEYS, "atmosphere")
    lapse = _value(t, "lapse_rate", Float64, "atmosphere")
    isfinite(lapse) || _fail("atmosphere.lapse_rate", "must be finite, got $lapse")
    return Atmosphere(;
        reference_speed = _nonnegative(
            _value(t, "reference_speed", Float64, "atmosphere"),
            "atmosphere.reference_speed",
        ),
        temperature = _positive(
            _value(t, "temperature", Float64, "atmosphere"),
            "atmosphere.temperature",
        ),
        density = _positive(
            _value(t, "density", Float64, "atmosphere"),
            "atmosphere.density",
        ),
        lapse_rate = lapse,
        specific_heat = _positive(
            _value(
                t,
                "specific_heat",
                Float64,
                "atmosphere";
                default = DRY_AIR_SPECIFIC_HEAT,
            ),
            "atmosphere.specific_heat",
        ),
        surface = _choice(
            _value(t, "surface", String, "atmosphere"),
            _SURFACE_CHOICES,
            "atmosphere.surface",
        ),
        roughness = _choice(
            _value(t, "roughness", String, "atmosphere"),
            _ROUGHNESS_CHOICES,
            "atmosphere.roughness",
        ),
    )
end

function _buildings_from(t::AbstractDict)
    _reject_unknown(t, _BUILDINGS_KEYS, "buildings")
    coefficient = _nonnegative(
        _value(
            t,
            "wake_coefficient",
            Float64,
            "buildings";
            default = DEFAULT_WAKE_COEFFICIENT,
        ),
        "buildings.wake_coefficient",
    )
    entries = get(t, "building", Any[])
    entries isa AbstractVector ||
        _fail("buildings.building", "expected an array of tables, got $(typeof(entries))")
    buildings = Building[]
    for (i, entry) in enumerate(entries)
        p = "buildings.building[$i]"
        entry isa AbstractDict || _fail(p, "expected a table, got $(typeof(entry))")
        _reject_unknown(entry, _BUILDING_KEYS, p)
        push!(
            buildings,
            Building(;
                east = _value(entry, "east", Float64, p),
                north = _value(entry, "north", Float64, p),
                height = _nonnegative(_value(entry, "height", Float64, p), "$p.height"),
                frontal_area = _nonnegative(
                    _value(entry, "frontal_area", Float64, p),
                    "$p.frontal_area",
                ),
            ),
        )
    end
    return BuildingEnvelope(buildings; wake_coefficient = coefficient)
end

function _nuclide_from(t::AbstractDict)
    _reject_unknown(t, _NUCLIDE_KEYS, "nuclide")
    low = _nonnegative(
        _value(t, "deposition_velocity_low", Float64, "nuclide"),
        "nuclide.deposition_velocity_low",
    )
    high = _nonnegative(
        _value(t, "deposition_velocity_high", Float64, "nuclide"),
        "nuclide.deposition_velocity_high",
    )
    low ≤ high || _fail(
        "nuclide.deposition_velocity_high",
        "must not be below nuclide.deposition_velocity_low ($low m/s), got $high",
    )
    return Nuclide(;
        name = _value(t, "name", String, "nuclide"),
        decay_constant = _nonnegative(
            _value(t, "decay_constant", Float64, "nuclide"),
            "nuclide.decay_constant",
        ),
        deposition_velocity = DepositionVelocity(low, high),
        washout_species = _choice(
            _value(t, "washout_species", String, "nuclide"; default = "tritium_iodine"),
            _SPECIES_CHOICES,
            "nuclide.washout_species",
        ),
    )
end

function _rose_from(t::AbstractDict)
    _reject_unknown(t, _ROSE_KEYS, "wind_rose")
    n = _value(t, "sectors", Int, "wind_rose"; default = 16)
    (n ≥ 4 && iseven(n)) ||
        _fail("wind_rose.sectors", "must be an even number of at least 4, got $n")
    grid = SectorGrid(n)

    convention = _choice(
        _value(t, "convention", String, "wind_rose"),
        _CONVENTION_CHOICES,
        "wind_rose.convention",
    )

    frequencies = _floatvector(t, "frequencies", "wind_rose")
    length(frequencies) == n || _fail(
        "wind_rose.frequencies",
        "must hold one value per sector ($n), got $(length(frequencies))",
    )

    stability = _floatvector(t, "stability", "wind_rose")
    nclasses = length(PASQUILL_CLASSES)
    length(stability) == nclasses || _fail(
        "wind_rose.stability",
        "must hold one value per Pasquill class ($nclasses), got $(length(stability))",
    )

    return WindRose(grid, frequencies, convention; stability)
end

function _floatvector(t::AbstractDict, key::AbstractString, path::AbstractString)
    haskey(t, key) || _fail(_join(path, key), "required key is missing")
    raw = t[key]
    raw isa AbstractVector ||
        _fail(_join(path, key), "expected an array of numbers, got $(typeof(raw))")
    out = Vector{Float64}(undef, length(raw))
    for (i, v) in enumerate(raw)
        v isa Real ||
            _fail("$(_join(path, key))[$i]", "expected a number, got $(typeof(v))")
        out[i] = Float64(v)
    end
    return out
end
