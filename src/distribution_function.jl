"""
$(SIGNATURES)

Mean of the initial condition in x.
"""
function mean_f0(v::Float64)::Float64
    return (1.0 / sqrt(2π)) * exp(-0.5 * v * v)
end

"""
$(SIGNATURES)
    
Initial condition for the Vlasov-equation (Penrose-stable equilibra with a perturbation).
"""
function f0(x::Float64, v::Float64)::Float64
    a = 0.01
    k = 1.0
    return (1.0 / sqrt(2π)) * exp(-0.5 * v * v) * (1 + a * cos(2π * k * x))
end

export compute_initial_condition

"""
$(SIGNATURES)
    
Compute the initial condition for rho, u, and total density.
"""
function compute_initial_condition(mesh)

    nx, ng = mesh.nx, mesh.ng
    rho = zeros(ng, nx + 1)
    u = zeros(ng, nx + 1)
    rho_tot = zeros(nx + 1)

    for j in 1:ng
        alpha = mesh.x[j]
        for i in 1:(nx+1)
            x_i = (i - 1) * mesh.dx
            rho[j, i] = f0(x_i, alpha) / mean_f0(alpha)
            u[j, i] = alpha
            rho_tot[i] += mesh.w[j] * rho[j, i] * mean_f0(alpha) * exp(alpha * alpha)
        end
    end

    return rho, u, rho_tot

end

