using DispersionRelations
using LinearAlgebra
using MultiStreamVlasovPoisson
using Plots
using .Threads
global k = 0.5          #Wave number
global u0 = 0.0         #Mean Velocity : 1) u0 = 0 Maxwellian, 2) u0 !=0 Two streams
global T = 1.0          #Temperature
global L  = 2π / k      #Size of the domain
global eps = 1.0        #Debye length

function main(hermite_quad)
    nx, xmin, xmax = 256, 0.0, L
    mesh_x = UniformMesh(xmin,xmax,nx)
    nv, vmin, vmax = 256, -6.0, 6.0 
    grid_v = UniformGrid(vmin, vmax, nv, T,u0)

    #Compute the initial condition
    rho, u, rho_tot = compute_initial_condition(mesh_x,grid_v,k,T,u0)
    phi = zeros(nx+1)
    poisson!(phi, mesh_x, rho_tot)

    #Initialize the streams
    u_at_step_n   = zeros(nx + 1,nv)    
    rho_at_step_n = zeros(nx + 1,nv)    

    #Set the CFL number and the final time
    dt =  0.1*mesh_x.dx
    tfinal = 50
    time = [0.0]

    #Array of physical quantities 
    elec_energy     = [compute_elec_energy(phi, mesh_x, eps)]
    mass            = [compute_total_mass(rho_tot,mesh_x)]
    momentum        = [compute_momentum(rho,u,mesh_x,grid_v)]
    total_energy    = [compute_elec_energy(phi, mesh_x, eps) + compute_kinetic_energy(rho,u,mesh_x,grid_v)]

    #Temporal loop
    n = 0
    while n * dt <= tfinal
	    iter = 0
        err=1e-10
	    maxiter=50
        norm_dx_u = 0.0
        #Compute the discrete L^{\infty} norm of the gradient of the velocities
        norm_dx_u = compute_norm_dx_u(mesh_x,grid_v,u)
        threshold = 1.0/(n*dt)
        if(norm_dx_u > threshold )
            println("Remapping f at time = $(n*dt)")
            rho, u = remap_f_on_uniform_grid(mesh_x,grid_v,rho,u)
        end
        copyto!(rho_at_step_n, rho)
        copyto!(u_at_step_n, u)	
        old_u = zeros(nx + 1,nv)
       
       #Fixed point lopp to solve the non linear MultiStream pressureless Euler-Poisson system
        while norm(u .- old_u, Inf) / norm(u, Inf) > err && iter < maxiter
            #Update rho : the streams are packed into groups that are solved on each threads
            @threads for j in 1:nv	 
                update_rho!(mesh_x, view(rho, :, j), view(u, :, j), view(rho_at_step_n, :, j), dt)
            end
            #Assemble rho
	        compute_rho_total!(rho_tot,grid_v,rho)
            #Solve Poisson
	        poisson!(phi, mesh_x, rho_tot)	  
            copyto!(old_u, u)
            #Update u : the streams are packed into groups that are solved on each threads
	        @threads for j in 1:nv
	  	        update_u!(mesh_x, view(rho, :, j), view(u, :, j), phi,
		        view(rho_at_step_n, :, j), view(u_at_step_n, :, j), dt, maxiter)
            end
         iter += 1
	    end
        compute_rho_total!(rho_tot, grid_v, rho)
    
        push!(elec_energy, compute_elec_energy(phi, mesh_x, eps))
        push!(mass,compute_total_mass(rho_tot,mesh_x))
        push!(momentum,compute_momentum(rho,u,mesh_x,grid_v))
        push!(total_energy,compute_elec_energy(phi, mesh_x, eps) + compute_kinetic_energy(rho,u,mesh_x,grid_v))
        n += 1
        push!(time, n * dt)
        println("iteration: $n , time = $(n * dt), elec energy = $(last(elec_energy)), mass = $(last(mass)),   ||dxU|| = $norm_dx_u")  
    end

    #Plot the final distribution function
    f_on_grid =  interpolate_f_on_grid(mesh_x,grid_v,rho,u)
    X = []
    Y = []
    Z = []
    for i in 1:nx+1
        for j in 1:nv
            X = push!(X,mesh_x.x[i])
            Y = push!(Y,grid_v.v[j])
            Z = push!(Z,f_on_grid[i,j])
        end
    end
    plot_f = plot(X,Y,Z,st = [:surface],camera = (0,90),xlabel = "x", ylabel="v")

    return time, elec_energy, mass, momentum, total_energy, grid_v, mesh_x, u, rho_tot, plot_f

end
@time time, elec_energy, mass, momentum, total_energy, grid_v, mesh_x, u, rho_tot, plot_f = main(true)
plot(time,elec_energy,yaxis = :log)

#plot(mesh_x.x,u[:,grid_v.nv/2-10:grid_v.nv/2+10], legend = false)

#p_1 = plot(time, elec_energy, yaxis = :log, label = "UniformGrid")
#p_2 = plot(X,Y,Z,st = [:surface], label="f")
#plot!(p_2,camera = (0,90))
# @time t, elec_energy = main(true)
#plot!(t, elec_energy, yaxis = :log, label = "hermite")

# line, ω, = fit_complex_frequency(t, elec_energy, use_peaks = 1:2)
# plot!(time, line; yaxis = :log)
# title!("ω = $(imag(ω))")


#for k=0.4, we have from Eric's book (a=0.001)
#E_k(t)=0.002.*0.424666.*exp.(-0.0661.*t).*abs.(cos.(1.285.*t.−0.3357725)))
#plot(t,log.(elec_energy))
#plot!(t,log.(0.002.*0.42466.*exp.(-0.0661.*t).*abs.(cos.(1.285.*t.-0.33577))))
#plot!(t,log.(0.0075.*0.42466.*exp.(-0.0661.*t).*abs.(cos.(1.285.*t.-0.33577))))   for a=0.001

#for k=0.5, we have from Eric's book
#E_k(t)=4.*a.*0.3677.*exp.(-0.1533.*t).*sqrt.(2pi).*abs.(cos(1.4156.*t-0.536245))
#plot(t,log.(elec_energy))
#plot!(t,log.(0.01.*exp.(-0.1533.*t).*sqrt(2*pi).*abs.(cos.(1.4156.*t.-0.536245))))           for a=0.01
#plot!(t,log.(0.0025.*0.3677.*exp.(-0.1533.*t).*sqrt(2*pi).*abs.(cos.(1.4156.*t.-0.536245))))  for a=0.001