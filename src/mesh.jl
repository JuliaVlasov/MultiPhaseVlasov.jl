using FastGaussQuadrature

abstract type AbstractMesh end

export GaussHermiteMesh

struct GaussHermiteMesh <: AbstractMesh

    nx::Int
    dx::Float64
    ng::Int
    x::Vector{Float64}
    w::Vector{Float64}

    function GaussHermiteMesh(nx, ng)

        dx = 1.0 / (nx + 1)
        x, w = gausshermite(ng)

        return new(nx, dx, ng, x, w)

    end

end

export UniformMesh

mutable struct UniformMesh <: AbstractMesh

    nx::Int
    dx::Float64
    ng::Int
    vmin::Float64
    vmax::Float64
    sf0::Float64

    function UniformMesh(nx, vmin, vmax, ng)

        dx = 1.0 / (nx + 1)
        sf0 = 0.0
        return new(nx, dx, ng, vmin, vmax, sf0)

    end

end
