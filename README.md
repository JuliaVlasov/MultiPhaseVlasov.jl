# Code VP Multistream

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://mingus.gitlabpages.inria.fr/code-vp-multistream/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://mingus.gitlabpages.inria.fr/code-vp-multistream/stable

[GHA-img]: https://gitlab.inria.fr/mingus/code-vp-multistream/badges/main/pipeline.svg
[GHA-url]: https://gitlab.inria.fr/mingus/code-vp-multistream/commits/main

## Sources

```
git clone git@gitlab.inria.fr:mingus/code-vp-multistream.git
cd code-vp-multistream
```

## Compilation C++

```
g++ -O3 -I. *.cpp -o main
```

## Run C++

```
./main
```

## Julia code

```
julia
julia> import Pkg
julia> Pkg.update()
julia> Pkg.add("Plots")
julia> Pkg.activate(pwd())
julia> Pkg.instantiate()
julia> include("main.jl")
```

