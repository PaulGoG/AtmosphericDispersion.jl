@testset "Configuration" begin
    reference = joinpath(@__DIR__, "..", "config", "reference.toml")

    @testset "the reference configuration loads" begin
        @test isfile(reference)
        c = load_configuration(reference)
        @test c isa RunConfiguration
        @test c.site.source.height == 50.3
        @test c.site.atmosphere.surface === SURFACE_AGRICULTURAL
        @test c.site.atmosphere.roughness === ROUGHNESS_PASTURE
        @test c.nuclide.name == "HTO"
        @test c.nuclide.decay_constant == TRITIUM_DECAY_CONSTANT
        @test nsectors(grid(c.rose)) == 16
        @test c.activity == 6.8e14
        @test c.washout.precipitation === PRECIPITATION_RAIN
        @test c.washout.rate == 1.0
        @test c.washout.duration == 0.0
        @test c.washout.model === WASHOUT_NORMATIVE
        @test c.spacing < c.extent
        # It drives the solver.
        @test dilution_long_term(0.0, -c.extent, c.site, c.rose; nuclide = c.nuclide) > 0
        @test_throws ArgumentError load_configuration(reference * ".missing")
    end

    # Mutate a copy of the reference table and require the named failure.
    base() = TOML.parsefile(reference)
    function withkey(f)
        t = base()
        f(t)
        return t
    end

    @testset "missing keys are named" begin
        for (table, key) in (
            ("source", "height"),
            ("source", "diameter"),
            ("atmosphere", "temperature"),
            ("atmosphere", "surface"),
            ("nuclide", "decay_constant"),
            ("release", "activity"),
            ("grid", "extent"),
            ("wind_rose", "convention")
        )
            t = withkey(t -> delete!(t[table], key))
            err = try
                configuration_from(t)
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "$table.$key"
            @test occursin("missing", err.message)
        end
        # A missing table is named too.
        t = withkey(t -> delete!(t, "release"))
        err = try
            configuration_from(t)
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "release"
    end

    @testset "enumerated choices list their options" begin
        for (table, key, options) in (
            ("atmosphere", "surface", "agricultural"),
            ("atmosphere", "roughness", "pasture"),
            ("wind_rose", "convention", "blowing_from"),
            ("precipitation", "type", "rain")
        )
            t = withkey(t -> (t[table][key] = "nonsense"))
            err = try
                configuration_from(t)
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "$table.$key"
            @test occursin("must be one of", err.message)
            @test occursin(options, err.message)
        end
    end

    @testset "numerical bounds" begin
        cases = [
            (t -> (t["source"]["height"] = 0.0), "source.height"),
            (t -> (t["source"]["diameter"] = -1.0), "source.diameter"),
            (t -> (t["source"]["exit_velocity"] = -1.0), "source.exit_velocity"),
            (t -> (t["atmosphere"]["reference_speed"] = 0.0), "atmosphere.reference_speed"),
            (t -> (t["atmosphere"]["temperature"] = 0.0), "atmosphere.temperature"),
            (t -> (t["atmosphere"]["density"] = -1.0), "atmosphere.density"),
            (t -> (t["release"]["activity"] = -1.0), "release.activity"),
            (t -> (t["release"]["duration"] = 0.0), "release.duration"),
            (t -> (t["grid"]["extent"] = -1.0), "grid.extent"),
            (t -> (t["grid"]["spacing"] = 0.0), "grid.spacing"),
            (
                t -> (t["buildings"]["wake_coefficient"] = -1.0),
                "buildings.wake_coefficient",
            ),
            (
                t -> (t["buildings"]["broadening_limit"] = 0.5),
                "buildings.broadening_limit",
            ),
            (
                t -> (t["precipitation"]["washout_duration"] = -1.0),
                "precipitation.washout_duration",
            ),
        ]
        for (mutate, path) in cases
            err = try
                configuration_from(withkey(mutate))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == path
        end
    end

    @testset "cross-key constraints" begin
        # The receptor spacing must fit inside the grid.
        err = try
            configuration_from(withkey(t -> (t["grid"]["spacing"] = t["grid"]["extent"])))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "grid.spacing"

        # The high deposition velocity cannot sit below the low one.
        err = try
            configuration_from(withkey(t -> (t["nuclide"]["deposition_velocity_high"] = 1e-6)))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "nuclide.deposition_velocity_high"

        # The washout table is defined only at tabulated intensities.
        err = try
            configuration_from(withkey(t -> (t["precipitation"]["rate"] = 2.0)))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "precipitation.rate"
    end

    @testset "the wind rose is checked on the way in" begin
        # One frequency per sector.
        err = try
            configuration_from(withkey(t -> (t["wind_rose"]["frequencies"] = fill(0.125, 8))))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "wind_rose.frequencies"

        # One stability fraction per Pasquill class.
        err = try
            configuration_from(withkey(t -> (t["wind_rose"]["stability"] = fill(0.25, 4))))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "wind_rose.stability"

        # Sector count must be even and at least four.
        for n in (3, 15, 2)
            err = try
                configuration_from(withkey(t -> (t["wind_rose"]["sectors"] = n)))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == "wind_rose.sectors"
        end

        # Raw counts rather than a distribution reach the WindRose constructor.
        @test_throws ArgumentError configuration_from(
            withkey(t -> (t["wind_rose"]["frequencies"] = fill(10.0, 16))),
        )
    end

    @testset "the declared convention reaches the field" begin
        from = configuration_from(
            withkey(t -> (t["wind_rose"]["convention"] = "blowing_from")),
        )
        toward = configuration_from(
            withkey(t -> (t["wind_rose"]["convention"] = "blowing_toward")),
        )
        g = grid(from.rose)
        for k in 1:16
            @test frequency_toward(from.rose, k) ==
                  frequency_toward(toward.rose, opposite(g, k))
        end
        χ(c, k) = dilution_long_term(
            5000sin(sector_bearing(g, k)),
            5000cos(sector_bearing(g, k)),
            c.site,
            c.rose,
        )
        @test argmax([χ(from, k) for k in 1:16]) ==
              opposite(g, argmax([χ(toward, k) for k in 1:16]))
    end

    @testset "types are enforced" begin
        err = try
            configuration_from(withkey(t -> (t["source"]["height"] = "tall")))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "source.height"
        @test occursin("expected", err.message)
        # An unadorned TOML integer is accepted where a float is wanted.
        @test configuration_from(withkey(t -> (t["source"]["height"] = 50))).site.source.height ==
              50.0
    end

    @testset "unknown keys are rejected, by their dotted path" begin
        # A key the loader does not read is indistinguishable from a
        # misspelt one, which would otherwise take its default in silence.
        cases = (
            (t -> t["sorce"] = Dict("height" => 1.0), "sorce"),
            (t -> t["source"]["hieght"] = 1.0, "source.hieght"),
            (t -> t["atmosphere"]["rougness"] = "pasture", "atmosphere.rougness"),
            (t -> t["model"]["mixing_layr"] = "unbounded", "model.mixing_layr"),
            (t -> t["buildings"]["wake_coeficient"] = 1.0, "buildings.wake_coeficient"),
            (t -> t["nuclide"]["half_life"] = 12.3, "nuclide.half_life"),
            (t -> t["wind_rose"]["sector"] = 16, "wind_rose.sector"),
            (t -> t["release"]["activty"] = 1.0, "release.activty"),
            (t -> t["precipitation"]["kind"] = "rain", "precipitation.kind"),
            (t -> t["grid"]["extend"] = 1.0, "grid.extend"),
            (
                t -> t["buildings"]["building"] = [
                    Dict(
                    "east" => 1.0,
                    "north" => 0.0,
                    "height" => 5.0,
                    "frontal_area" => 10.0,
                    "width" => 3.0,
                ),
                ],
                "buildings.building[1].width",
            ),
        )
        for (mutate, path) in cases
            err = try
                configuration_from(withkey(mutate))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == path
            @test occursin("unknown key", err.message)
        end
        # An optional table of the wrong type is named rather than crashing.
        for table in ("buildings", "precipitation", "model")
            err = try
                configuration_from(withkey(t -> t[table] = 3))
                nothing
            catch e
                e
            end
            @test err isa ConfigurationError
            @test err.path == table
        end
    end

    @testset "the washout species is read from the nuclide table" begin
        @test load_configuration(reference).nuclide.washout_species ===
              WASHOUT_TRITIUM_IODINE
        c = configuration_from(withkey(t -> t["nuclide"]["washout_species"] = "other"))
        @test c.nuclide.washout_species === WASHOUT_OTHER_NUCLIDES
        err = try
            configuration_from(withkey(t -> t["nuclide"]["washout_species"] = "caesium"))
            nothing
        catch e
            e
        end
        @test err isa ConfigurationError
        @test err.path == "nuclide.washout_species"
    end
end

@testset "Configuration: dispersion scheme, stable wind and extrapolation" begin
    reference = joinpath(@__DIR__, "..", "config", "reference.toml")
    withkey(f) = (t = TOML.parsefile(reference); f(t); t)
    failure(t) =
        try
            configuration_from(t)
            nothing
        catch e
            e
        end

    c = configuration_from(withkey(t -> begin
        t["model"]["dispersion"] = "briggs_urban"
        t["model"]["stable_rise_wind"] = "release_height"
    end))
    @test dispersion_scheme(c.site) === DISPERSION_BRIGGS_URBAN
    @test c.site.stable_rise_wind === WIND_AT_RELEASE_HEIGHT
    default = load_configuration(reference)
    @test dispersion_scheme(default.site) === DISPERSION_HOSKER
    @test default.site.stable_rise_wind === WIND_MEAN_OVER_RISE
    for (key, option) in (
        ("dispersion", "briggs_open_country"),
        ("stable_rise_wind", "mean_over_rise"),
        ("extrapolation", "refuse")
    )
        err = failure(withkey(t -> (t["model"][key] = "nonsense")))
        @test err isa ConfigurationError
        @test err.path == "model.$key"
        @test occursin(option, err.message)
    end

    # The reference grid reaches 20 km, twice the band of every Briggs-based
    # scheme, and the file sets "allow" to accept that.
    @test_logs configuration_from(withkey(identity))
    @test_logs (:warn, r"validity range") configuration_from(
        withkey(t -> (t["model"]["extrapolation"] = "warn")),
    )
    err = failure(withkey(t -> (t["model"]["extrapolation"] = "refuse")))
    @test err isa ConfigurationError
    @test err.path == "grid.extent"
    @test occursin("20000.0 m lies beyond", err.message)
    # Inside the band every policy is silent.
    @test_logs configuration_from(withkey(t -> begin
        t["model"]["extrapolation"] = "refuse"
        t["grid"]["extent"] = 10_000.0
        t["grid"]["spacing"] = 100.0
    end))
    # The nearest receptor is checked too, on its own key.
    err = failure(withkey(t -> begin
        t["model"]["extrapolation"] = "refuse"
        t["grid"]["spacing"] = 50.0
        t["grid"]["extent"] = 5000.0
    end))
    @test err isa ConfigurationError
    @test err.path == "grid.spacing"
    # Eimutis–Konicek reaches 100 km, so the reference grid lies inside its band.
    @test_logs configuration_from(withkey(t -> begin
        t["model"]["extrapolation"] = "refuse"
        t["model"]["dispersion"] = "eimutis_konicek"
    end))
end

@testset "Configuration: precipitation" begin
    reference = joinpath(@__DIR__, "..", "config", "reference.toml")
    withkey(f) = (t = TOML.parsefile(reference); f(t); t)
    failure(t) =
        try
            configuration_from(t)
            nothing
        catch e
            e
        end
    # The loader passes the release duration to the event: an hour of rain in
    # the year-long reference release is that share of it. So small a share
    # leaves the airborne fraction nearly the same whenever the hour falls.
    c = @test_logs configuration_from(
        withkey(t -> (t["precipitation"]["washout_duration"] = 3600.0)),)
    @test precipitation_fraction(c.washout) ≈ 3600 / 31_557_600
    # Durations, activity and extent must be finite.
    for (table, key) in (("precipitation", "washout_duration"), ("release", "duration"),
        ("release", "activity"), ("grid", "extent"))
        err = failure(withkey(t -> (t[table][key] = Inf)))
        @test err isa ConfigurationError
        @test err.path == "$table.$key"
        @test occursin("finite", err.message)
    end
    # The loader warns when a precipitation beginning with the release would
    # leave more than 10 % more airborne at the edge of the grid, in the
    # slowest class (4254 s to 20 km), than the washout form. None, or a long
    # one, is not reported, and neither is weak scavenging: an hour of rain on
    # HTO through a one-hour release (Λ = 10⁻⁵ s⁻¹) changes it by 2.5 %.
    @test_logs configuration_from(withkey(identity))
    @test_logs configuration_from(
        withkey(t -> (t["precipitation"]["washout_duration"] = 86_400.0)),)
    @test_logs configuration_from(withkey(t -> begin
        t["precipitation"]["washout_duration"] = 3600.0
        t["release"]["duration"] = 3600.0
    end))
    # Where Λ x/u is large outlasting the travel time is not enough: 10 h of
    # snow at 5 mm/h on a species of the other row (Λ = 10⁻³ s⁻¹) through a
    # 10 h release leave 0.0399 airborne at 20 km against 0.0142.
    @test_logs (:warn, r"0\.0399.*0\.0142") configuration_from(withkey(t -> begin
        t["precipitation"]["washout_duration"] = 36_000.0
        t["release"]["duration"] = 36_000.0
        t["precipitation"]["type"] = "snow"
        t["precipitation"]["rate"] = 5.0
        t["nuclide"]["washout_species"] = "other"
    end))
    # Precipitation that outlasts the release leaves the most airborne when it
    # ends with the release: 8000 s of that snow over a 3000 s release, which
    # the form puts entirely in it.
    snow(t) = (t["precipitation"]["type"] = "snow"; t["precipitation"]["rate"] = 5.0;
        t["nuclide"]["washout_species"] = "other")
    @test_logs (:warn, r"washout_duration") configuration_from(withkey(t -> begin
        snow(t)
        t["precipitation"]["washout_duration"] = 8000.0
        t["release"]["duration"] = 3000.0
    end))
    # The warning is formed for any finite input, also where the form underflows.
    for (extent, t_w, T_R) in ((1e20, 1e19, 1e19), (1e20, 1e19, 1e18), (5e6, 1e6, 1e5))
        @test_logs (:warn, r"washout_duration") configuration_from(withkey(t -> begin
            snow(t)
            t["grid"]["extent"] = extent
            t["precipitation"]["washout_duration"] = t_w
            t["release"]["duration"] = T_R
        end))
    end
end
