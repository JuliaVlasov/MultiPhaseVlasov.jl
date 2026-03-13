# Vlasov-Poisson with Multistream numerical method

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliavlasov.github.io/MultiStreamVlasovPoisson.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliavlasov.github.io/MultiStreamVlasovPoisson.jl/stable

```
git clone https://github.com/JuliaVlasov/MultiStreamVlasovPoisson.jl.git
cd MultiStreamVlasovPoisson.jl
```

```
julia
julia> import Pkg
julia> Pkg.update()
julia> Pkg.add("Plots")
julia> Pkg.activate(pwd())
julia> Pkg.instantiate()
julia> include("main.jl")
```

