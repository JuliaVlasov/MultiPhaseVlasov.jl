export compute_rho_total!

"""
$(SIGNATURES)
    
Compute the total density from the density matrix.
"""
function compute_rho_total!(rho_tot::Vector{Float64}, mesh::GaussHermiteMesh, rho)

    fill!(rho_tot, 0.0)

    for i in eachindex(rho_tot)
        for j in axes(rho, 2)
            alpha = mesh.x[j]
            rho_tot[i] += mesh.w[j] * rho[i, j] * mean_f0(alpha) * exp(alpha * alpha)
        end
    end

    return
end
