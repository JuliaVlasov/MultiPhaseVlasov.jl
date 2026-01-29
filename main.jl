using LinearAlgebra
using MultiStreamVlasovPoisson
using Plots
using .Threads

function main(hermite_quad)

    eps = 1.0
    nx = 100
    vmin, vmax = -4.0, 4.0
    ng = 100

    if hermite_quad 
        mesh = GaussHermiteMesh(nx, ng)
    else
        mesh = UniformMesh(nx, vmin, vmax, ng)
    end

    rho, u, rho_tot = compute_initial_condition(mesh)

    poisson = NonLinearPoissonSolver(eps, nx)

    phi = -log.(rho_tot)

    dt = mesh.dx
    tfinal = 100 * dt  # Final time
    time = [0.0]

    elec_energy = [compute_elec_energy(phi, mesh, eps)]

    n = 0
    while n * dt <= tfinal

        # Update phi
        solve!(phi, poisson, rho_tot)
        for j in 1:ng
            update_single_fluid_solution!(mesh, poisson, view(rho, :, j), view(u, :, j), rho_tot, phi, dt)
        end
        compute_rho_total!(rho_tot, mesh, rho)
        push!(elec_energy, compute_elec_energy(phi, mesh, eps))
        n += 1
        push!(time, n * dt)
        println("iteration: $n , time = $(n * dt), elec energy = $(last(elec_energy))")

    end

    return time, elec_energy

end

@time time, elec_energy = main(false)
plot(time, elec_energy, yscale = :log10, label = "uniform")

@time time, elec_energy = main(true)
plot!(time, elec_energy, yscale = :log10, label = "hermite")


