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
            ("wind_rose", "convention"),
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
            ("precipitation", "type", "rain"),
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
        for k = 1:16
            @test frequency_toward(from.rose, k) ==
                  frequency_toward(toward.rose, opposite(g, k))
        end
        χ(c, k) = dilution_long_term(
            5000sin(sector_bearing(g, k)),
            5000cos(sector_bearing(g, k)),
            c.site,
            c.rose,
        )
        @test argmax([χ(from, k) for k = 1:16]) == opposite(g, argmax([χ(toward, k) for k = 1:16]))
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
                t ->
                    t["buildings"]["building"] = [
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
