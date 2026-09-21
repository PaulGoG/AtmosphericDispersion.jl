using Pkg
Pkg.activate(@__DIR__; io = devnull)
# [sources] is read from Julia 1.11 on; before that the package is developed by path.
VERSION < v"1.11" && Pkg.develop(; path = dirname(@__DIR__), io = devnull)
Pkg.instantiate(; io = devnull)
