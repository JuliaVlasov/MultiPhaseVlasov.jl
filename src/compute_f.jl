export interpolate_f_on_grid
export remap_f_on_uniform_grid
using .Threads

function interpolate_f_on_grid(mesh_x::AbstractMesh, grid_v::UniformGrid,rho::Matrix{Float64},u::Matrix{Float64})
    nx = mesh_x.nx
    nv = grid_v.nv
    dv = grid_v.dv
    rf = zeros(nx+1,nv)
    for i in 1: (nx+1)
        for j in 1:nv
            for l in 1:nv
                v_j = grid_v.v[j]
                rf[i,j] += grid_v.w[l] * rho[i,l] * Spline(v_j-u[i,l],dv) ##Naive approchae we do not localize v_j-u_{i,l} to optimize later
            end
        end
    end
    return rf
end
#This works only for uniform grid for the moment
function remap_f_on_uniform_grid(mesh_x::AbstractMesh, grid_v::UniformGrid, rho::Matrix{Float64}, u::Matrix{Float64})
    nv = grid_v.nv
    dv = grid_v.dv
    nx = mesh_x.nx
    new_SF = 0.0
    new_weights = zeros(nv)
    new_rho = zeros(nx+1, nv)
    new_u   = zeros(nx+1, nv)
    #Reinitiaize rho and u on uniform grid
    @threads for l in 1:nv
        for i in 1:(nx+1)
            new_u[i,l] = grid_v.v[l]
            new_rho[i,l] = evaluate_f_on(i,l,grid_v,rho,u)/ evaluate_mean_f_on(l,mesh_x,grid_v,rho,u)
        end
    end
    #Compute the new weights
    @threads for l in 1:nv
        new_weights[l] = evaluate_mean_f_on(l,mesh_x,grid_v,rho,u) * dv
        new_SF+= new_weights[l]
    end
    #Update the weights
    @threads for l in 1:nv
        new_weights[l] *=(1.0/new_SF)
        grid_v.w[l] = new_weights[l]
    end
    return new_rho, new_u
end
#Evaluate f on (x_i,v_j)
function evaluate_f_on(i::Int,j::Int,grid_v::UniformGrid, rho::Matrix{Float64}, u::Matrix{Float64})
    nv = grid_v.nv
    dv = grid_v.dv
    v_j = grid_v.v[j]
    f = 0.0
    for l in 1:nv
        f+= grid_v.w[l] * rho[i,l] * Spline(v_j - u[i,l],dv)
    end
    return f
end
#Compute mean_f at v_j on the torus of size L = xmax - xmin
function evaluate_mean_f_on(j::Int,mesh_x::AbstractMesh,grid_v::UniformGrid,rho::Matrix{Float64},u::Matrix{Float64})
    mean_f = 0.0
    nv = grid_v.nv
    dv = grid_v.dv
    nx = mesh_x.nx
    dx = mesh_x.dx
    v_j = grid_v.v[j]
    for  l in 1:nv
        for i in 1:(nx+1)
            mean_f+= (grid_v.w[l] * rho[i,l] * Spline(v_j-u[i,l],dv) * dx) /(mesh_x.x[nx+1]- mesh_x.x[1])
        end
    end
    return mean_f
end



#Spline of order 2 supp S_2 = [-1.5 h , 1.5 h]
function Spline(x::Float64,h::Float64)
    if( 0.5 * h < abs(x) < 1.5 * h)
        S = (0.5 * (1.5-abs(x/h))*(1.5-abs(x/h)))/h
    elseif ( abs(x) < 0.5 * h)
        S = ((3.0/4.0) -(x/h)*(x/h))/h
    else
        S = 0
    end
    return S
end
