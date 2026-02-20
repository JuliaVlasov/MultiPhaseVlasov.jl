export mean_f0
export f0
"""
# Here is defined the initial condition:
1) Maxwellian case is when (u_0 = 0)
2) Two stream case is when (u_0 != 0)
$(SIGNATURES)
Mean of the initial condition in x // The perturbation in x must be of zero mean
"""
function mean_f0(v::Float64, T::Float64,u0::Float64)::Float64
    return (1.0 / sqrt(2π*T)) *   ( 0.5 * exp(-0.5 * (v-u0) * (v-u0)/T) + 0.5* exp(-0.5 * (v+u0) * (v+u0)/T) ) 
end

"""
$(SIGNATURES)
Initial condition for the Vlasov-equation (Penrose-stable equilibra with a perturbation).
"""
function f0(x::Float64, v::Float64, k::Float64, T::Float64, u0::Float64)::Float64
    a = 0.01
    return (1.0 / sqrt(2π*T)) *   ( 0.5 * exp(-0.5 * (v-u0) * (v-u0)/T) + 0.5* exp(-0.5 * (v+u0) * (v+u0)/T) ) * (1 + a * cos(k * x))
end



