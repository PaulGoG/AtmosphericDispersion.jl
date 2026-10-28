# The worked run of the manual and a configuration as scripts/run.jl reads one,
# compiled into the package image so that the first run of a session does not
# spend seconds in compilation.
@setup_workload begin
    configuration = """
    [source]
    height = 50.3
    diameter = 2.33
    exit_velocity = 10.0
    exit_density = 0.6
    exit_temperature = 324.0
    [atmosphere]
    reference_speed = 4.0
    temperature = 287.0
    density = 1.2
    lapse_rate = 0.02
    specific_heat = 1005.0
    surface = "agricultural"
    roughness = "pasture"
    [model]
    plume_rise = "briggs"
    stable_rise_wind = "mean_over_rise"
    dispersion = "hosker"
    extrapolation = "allow"
    resuspension = "iaea_ss57"
    washout = "normative"
    [mixing_layer]
    scheme = "tabulated"
    above_lid = "rise_inhibited"
    [buildings]
    wake_coefficient = 1.0
    [nuclide]
    name = "HTO"
    decay_constant = 1.784e-9
    deposition_velocity_low = 0.004
    deposition_velocity_high = 0.008
    washout_species = "tritium_iodine"
    [wind_rose]
    sectors = 4
    convention = "blowing_from"
    frequencies = [0.25, 0.25, 0.25, 0.25]
    stability = [0.1, 0.1, 0.1, 0.4, 0.15, 0.15]
    [release]
    activity = 1.0e14
    duration = 3.15576e7
    [precipitation]
    type = "rain"
    rate = 1.0
    washout_duration = 0.0
    [grid]
    extent = 5000.0
    spacing = 1000.0
    """
    @compile_workload begin
        config = configuration_from(TOML.parse(configuration))
        (; nuclide, washout) = config
        dilution_long_term(0.0, -config.extent, config.site, config.rose; nuclide, washout)
        source = StackSource(;
            height = 50.3, diameter = 2.33, exit_velocity = 10.0,
            exit_density = 0.6, exit_temperature = 324.0,
        )
        atmosphere = Atmosphere(;
            reference_speed = 4.0, temperature = 287.0, density = 1.2,
            lapse_rate = 0.02, surface = SURFACE_AGRICULTURAL,
            roughness = ROUGHNESS_PASTURE,
        )
        site = Site(; source, atmosphere)
        stability = fill(1 / 6, 6)
        rose = WindRose(SectorGrid(16), fill(1 / 16, 16), BlowingFrom(); stability)
        effective_height(1000.0, site, PASQUILL_D)
        χQ = dilution_long_term(0.0, -5000.0, site, rose; nuclide = TRITIATED_WATER)
        deposition = dry_deposition(χQ, TRITIATED_WATER)
        resuspended_concentration(deposition, 1.0, RESUSPENSION_IAEA_SS57)
    end
end
