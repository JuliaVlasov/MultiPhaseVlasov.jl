using MultiStreamVlasovPoisson
using Plots

function main()

    eps = 1.0
    nx = 100
    vmin, vmax = -4.0, 4.0
    ng = 100

    mesh = GaussHermiteMesh(eps, nx, ng)

    rho, u, rho_tot = compute_initial_condition(mesh)

    poisson = NonLinearPoissonSolver(eps, nx)

    phi = -log.(rho_tot)

    solve!(phi, poisson, rho_tot)

    dt = mesh.dx
    tfinal = 1000 * dt  # Final time
    time = [0.0]
    
    elec_energy = [compute_elec_energy(phi, mesh)]
    
    n = 0
    while n * dt <= tfinal

        # for j in 1:ng
        #     update_single_fluid_sol(rho[j, :], u[j, :], rho_tot, phi, dt)
        # end
            
        compute_rho_total!(rho_tot, mesh, rho)
        push!(elec_energy, compute_elec_energy(phi, mesh))
        n += 1
        push!(time, n * dt)
        println("iteration: $n and time = $(n * dt)")
            
    end

    return time, elec_energy

end
    

@time time, elec_energy =  main()

plot(time, elec_energy)
