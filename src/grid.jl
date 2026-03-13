using Sobol
using FastGaussQuadrature
abstract type AbstractGrid end
export UniformGrid

"""
$(SIGNATURES) 
Construct the atomic measure approximation dmu_0^{n_v} 
Uses uniform mesh of the interval [vmin,vmax] with nv point with a rectangle formula for the weights.
"""

struct UniformGrid <: AbstractGrid
    nv::Int
    v::Vector{Float64}
    dv::Float64
    w::Vector{Float64}
    T::Float64

    function UniformGrid(vmin, vmax, nv,T,mesh_x,test_case)
        
        sf0 = 0.0
        dv = Float64
        dv = (vmax-vmin)/(nv)
        v = LinRange(vmin,vmax,nv)
        w = zeros(nv)
        for i in 1:nv
            sf0+= mean_f0(v[i],T,mesh_x,test_case) * dv  
        end
        for i in 1:nv
            w[i] = mean_f0(v[i],T,mesh_x,test_case) * dv / (sf0)
        end
        new(nv, v, dv, w, T)
    end

end

export  GaussHermiteGrid

"""
$(SIGNATURES) 
Construct the atomic measure approximation dmu_0^{n_v} 
Uses the Gauss Hermite quadrature points 

"""
struct GaussHermiteGrid <: AbstractGrid
    nv::Int
    v::Vector{Float64}
    dv::Float64
    w::Vector{Float64}
    T::Float64
    u0::Float64

    function GaussHermiteGrid(nv,T,mesh_x,test_case)
        v , w = gausshermite(nv) 
        dv = 1.0
        for i in 1:nv
            v[i] = sqrt(2*T) * v[i]
            w[i] = w[i] * mean_f0(v[i],T,mesh_x,test_case) * exp(v[i]*v[i]/(2*T)) * sqrt(2*T)
            if(i<= nv-1)
                dv = min(dv,abs(v[i+1]-v[i]))
            end
        end
        new(nv,v,dv,w,T,u0)
    end

end

#Function to sample a Maxwellian with mean and temperature given
 function sample_maxwellian(T::Float64, mean::Float64, nb_sample::Int64)

    v = Float64[]
    rand_1 = rand(1,nb_sample)
    rand_2 = rand(1,nb_sample)
    for i = 1:nb_sample
        r = sqrt(T) *sqrt(-2*log(rand_1[i])) + mean
        t = 2*π*rand_2[i]
        push!(v,r*cos(t))
    end
    return v
 end

 """
$(SIGNATURES) 
Construct the atomic measure approximation dmu_0^{n_v} 

Uses the Monte Carlo approach : sample a Maxwellian with given temperature and mean 

"""

export MonteCarloGrid
struct MonteCarloGrid <:AbstractGrid
    nv::Int
    v::Vector{Float64}
    w::Vector{Float64}
    u0::Float64
    function MonteCarloGrid(nv,T,u0,mesh_x,test_case)
        v = sample_maxwellian(T,u0,nv)
        v = sort(v)
        w = zeros(nv)
        sf0 = 0.0
        for i in 1:nv-1
            sf0+= mean_f0(v[i],T,mesh_x,test_case) * (v[i+1]-v[i])
        end
        for i in 1:nv-1
            w[i] = mean_f0(v[i],T,mesh_x,test_case)*(v[i+1]-v[i])/sf0
        end
        new(nv,v,w,u0)
    end
end

