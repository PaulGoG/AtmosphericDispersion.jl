```@meta
CurrentModule = AtmosphericDispersion
```

# Configuration

A run is described by a TOML file, never by editing source. The loader checks
presence, type, enumerated choice and numerical bound, and attaches every
failure to the key it failed on by its full dotted path.

```julia-repl
julia> load_configuration("broken.toml")
ERROR: ConfigurationError at `grid.spacing`: must be smaller than grid.extent (20000.0 m), got 20000.0
```

Validation runs at load and constructs the solver types immediately, so
anything those constructors reject in turn — frequencies that do not sum to
one, a negative stack height — fails at the same point rather than deep into a
run.

```@autodocs
Modules = [AtmosphericDispersion]
Pages = ["config.jl"]
```
