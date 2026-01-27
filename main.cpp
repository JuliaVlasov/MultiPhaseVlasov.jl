#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include "Matrix.h"
#include "quadrature.hpp"
using namespace std;


const double PI     = 4*atan(1.0);
const double  eps   =  1.0;
const int  f_eps    = floor(1./eps);
const int Nx         = 100; // Number of mesh points
const double dx     = 1./(Nx+1);
double vmin = -4;
double vmax = 4;
double SF0(0);
int nG       = 100; // The number of points in the Gaussian-Hermite quadrature in velocity
const quadrature quad = quadrature(nG); // Load the Gauss-Hermite quadrature
const bool hermite_quad = true; // Using hermite node in velocity otherwise use uniform grid in velocity not implemented yet !
/*
 This code solves the Vlasov-Poisson equation using a multi-stream formulation :
 1) Each stream verifies a pressure-less Euler-Poisson-Boltzmann equations on the torus.
 2) The numerical discretization of the fluid equations uses a staggered discretization where the mass flux is regularized. It was proven to be uncondtionally stable and AP.
 This is code is the property of : Mehdi Badsi, Daniel Han-Kwan and Nicolas Crouseilles.
 
 */
int main(int argc, char ** argv)
{
  //Prototypes
  
  int mod(int,int);
  void compute_initial_condition(Matrix &, Matrix &, Vector&);
  void compute_rho_total(Matrix &, Vector &);
  void solve_non_linear_Poisson(Vector &, Vector &, int = 1000);
  void update_single_fluid_sol(Vector &, Vector &, Vector &, Vector &, double, int  = 1000);
  void compute_elec_energy(Vector &, double &);
  void plot_sol(int, Matrix &, Matrix &, Vector &);
  double dt = dx;
  double T = 1000*dt; // Final Time
  double time(0);
  double elec_energy(0);
    
    
  if(hermite_quad == false){ nG = 200 ;} // Number of Points in Velocity in the uniform grid
    
    
  Matrix rho = Matrix(nG,Nx+1);
  Matrix u   = Matrix(nG,Nx+1);
  Vector rho_tot = Vector(Nx+1);
  Vector phi = Vector(Nx+1);
    
  //Compute the initial condition
    compute_initial_condition(rho,u,rho_tot);
    plot_sol(0,rho,u,phi);
    //Initialize Newton
    for(int i= 0; i <= Nx ; i++)
    {
        phi[i] = -log(rho_tot[i]);
    }
  
    // Initialize Poisson
    solve_non_linear_Poisson(rho_tot,phi);
    compute_elec_energy(phi,elec_energy);
    //Plot the diag
    string file_name = "sim/elec_energy.dat";
    ofstream file(file_name);
    file << 0 << "  " << elec_energy << endl;
    //Emod_0 = Emod;
   //Temporal Loop
    int n(0);
    while (n * dt <=T)
    {
        for(int j = 0 ; j < nG ; j++)
        {
            //cout << "rho[j] = " << rho.operator[](j) << " size of rho[j] == " << rho.operator[](j).getSize() << endl;
            update_single_fluid_sol(rho.operator[](j),u.operator[](j),rho_tot,phi,dt);
        }
        compute_rho_total(rho,rho_tot);
        plot_sol(n+1,rho,u,phi);
        compute_elec_energy(phi,elec_energy);
        cout << "Plot sol at iteration  :" << n+1 << " and time = " << n * dt << endl;
        file << n*dt << "  " << elec_energy << endl;
        n++;
    }
    file.close();
    
  return 0;
 
}
// The initial condition for the Vlasov-equation (Penrose-stable equilibra with a perturbation)
double f0(double x,double v)
{
    double a = 0.01;
    double k = 1.0;
    return (1./sqrt(2*PI)) * exp(-0.5 * v *v) * ( 1 + a * cos(2*PI*k*x));
}
double mean_f0(double v)
{
    return (1./sqrt(2*PI)) * exp(-0.5 * v * v);
}

void compute_initial_condition(Matrix & rho, Matrix & u, Vector & rho_tot)
{
    double x_i;
    double alpha; double w_j(0); double S(0);
    if(hermite_quad == true)
    {
        for(int j = 0 ; j < nG; j++)
        {
            alpha = quad._X[j];
            
            for(int i = 0 ; i <= Nx; i ++)
            {
                x_i = i * dx;
                rho[j][i] = f0(x_i,alpha)/mean_f0(alpha);
                u[j][i]  = alpha;
                rho_tot[i] += quad._W[j] * rho[j][i] * mean_f0(alpha) * exp(alpha*alpha);
            }
        }
    }
    else
    {
        for(int j = 0 ; j < nG; j++)
        {
            alpha = vmin + j*(vmax-vmin)/(nG-1);
            SF0+= mean_f0(alpha);
        }
        for(int j = 0 ; j < nG; j++)
        {
            alpha = vmin + j*(vmax-vmin)/(nG-1);
            for(int i = 0 ; i <= Nx; i ++)
            {
                x_i = i * dx;
                rho[j][i] = f0(x_i,alpha)/mean_f0(alpha);
                u[j][i]  = alpha;
                rho_tot[i] += rho[j][i] * mean_f0(alpha)/SF0;
            }
        }
    }
  
    
}
void compute_rho_total(Matrix & rho, Vector & rho_tot)
{
    double x_i;
    double alpha; double S(0);
    if(hermite_quad == true)
    {
        for(int i = 0 ; i <= Nx ; i++)
        {
            x_i = i * dx;
            for(int j = 0 ; j < nG ; j++)
            {
                alpha = quad._X[j];
                rho_tot[i] += quad._W[j] * rho[j][i] * mean_f0(alpha) * exp(alpha*alpha);
            }
        }
    }
    else
    {
            for(int j = 0 ; j < nG; j++)
            {
                alpha = vmin + j*(vmax-vmin)/(nG-1);
                for(int i = 0 ; i <= Nx; i ++)
                {
                    x_i = i * dx;
                    rho_tot[i] += rho[j][i] * mean_f0(alpha)/SF0;
                }
            }
    }
}

double max(double a, double b)
{
    if (a>= b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

// Function g : here it is some regularization of the (x)^{+}.
double g(double x)
{
    //return max(0,x);
    if( x < -dx)
    {
        return 0;
    }
    else if( x >= -dx && x <= dx)
    {
        return (x+dx)*(x+dx)/(4*dx);
    }
    else
    {
        return x;
    }
}
double dg(double x)
{
    double h = 1E-16;
    return (g(x+h)-g(x-h))/(2*h);
}
// Function G : (s,t,u) --> s g(u) - t(g(u)-u)
double G(double s, double t, double u)
{
    return s * g(u) - t* (g(u)-u);
}
double du_G(double s, double t, double u)
{
    return -(t-s) * dg(u);
}
double t_rho(double s,double t, double u)
{
    if( abs(u) < 1E-16)
    {
       return t - (t-s) * dg(0);
    }
   else
   {
       return (G(s,t,u)-G(s,t,0))/u;
    }
}

int mod(int a,int b)
{
  int r = a%b;
  if(r < 0) r+=b;
  return r;
}


void solve_non_linear_Poisson(Vector & rho_tot, Vector & phi, int itermax = 1000)
{
    //This function solves the quasi-linear Elliptic problem: epsilon^2 Delta(phi) + e(-phi) = rho on the Torus
    // We use a Newton Method : We set F(phi)_i = espilon^2 Delta(phi)_i +e(-phi_i) - rho_i, i = 0,...,N
    //cout << "           Call : solve_non_linear_Poisson." << endl;
    int iter(0); int ok;
    Vector delta = Vector(Nx+1); Vector P = Vector(Nx+1);
    Vector F = Vector(Nx+1); F[0] = 1;
    Matrix Jac      = Matrix(Nx+1,Nx+1);
    while(iter < itermax && F.linfty() >= 1E-7 )
    {
        //Compute  F : R^{N+1} ----> R^{N+1} given above and its Jacobian Matrix
        F[0]        = eps * eps * (phi[1]- 2*phi[0] + phi[Nx]) + exp(-phi[0]) - rho_tot[0];
        Jac[0][0]   = -2 * eps * eps - exp(-phi[0]);
        Jac[0][1]   = eps * eps;
        Jac[0][Nx] = eps * eps;
        for(int i = 1 ; i <= Nx-1;i++)
        {
            F[i]        = eps * eps * (phi[i+1] -2 * phi[i] + phi[i-1]) + exp(-phi[i]) - rho_tot[i];
            Jac[i][i-1] = eps * eps;
            Jac[i][i]   = - 2 * eps * eps - exp(-phi[i]);
            Jac[i][i+1] = eps * eps;
            
        }
        F[Nx]        = eps * eps * (phi[0] - 2 * phi[Nx] + phi[Nx-1]) + exp(-phi[Nx]) - rho_tot[Nx];
        Jac[Nx][0]   = eps * eps;
        Jac[Nx][Nx]   = - 2 * eps * eps - exp(-phi[Nx]);
        Jac[Nx][Nx-1] = eps * eps;
        
        // Solve Jac F Delta = S
        //Factorize LU
        for(int i = 0; i <= Nx ;i++)P[i] = i;
        ok = Gauss(Jac,P);
        if(ok != 0){delta = Solve(Jac,P,F);}
        // Update the sol
        phi = phi - delta;
        //new_rho = phi;
        //cout << "                   iter_poisson  = " << iter << "     | F(phi^{"<<iter<<"}) | = " << F.linfty() << endl;
        // Reinitialize Jac and P an delta
        for(int i =0 ; i <=Nx;i++)
        {
            P[i] = i; delta[i] = 0;
            for(int j = 0 ; j<= Nx ; j++)
            {
                Jac[i][j] = 0;
            }
        }
        iter++;
        
        
    }
    //cout << "           EndCall : solve_non_linear_Poisson." << endl;
    
}

void compute_elec_energy(Vector & phi, double & E)
{
    E = 0;
    int ir(0);
    for(int i = 0; i <= Nx; i++)
    {
        ir = mod(i+1,Nx+1);
        E         += 0.5 * eps*eps * dx* (phi[ir]-phi[i])*(phi[ir]-phi[i])/(dx*dx);
    }
}
//Update one single fluid sol
void update_single_fluid_sol(Vector & rho, Vector & u , Vector & rho_tot, Vector & phi, double dt, int itermax = 1)
{
    // This function solves the discrete pressure-less Euler-Poisson equation by a fixed point
    //cout << "   Call : solve_pressure-less_Euler-Poisson." << endl;
    int iter(0); int ok; int ir; int il; int irr; int ill; double Q_r(0); double Q_l(0);
    Vector rho_at_step_n = Vector(Nx+1); rho_at_step_n = rho;
    Vector u_at_step_n   = Vector(Nx+1); u_at_step_n = u;
    Vector new_rho  = Vector(Nx+1); Vector old_u = Vector(Nx+1);
    Vector new_phi = Vector(Nx+1); Vector q = Vector(Nx+1);
    Matrix M_upd_rho  = Matrix(Nx+1,Nx+1); Matrix M_upd_u = Matrix(Nx+1,Nx+1);  Vector P = Vector(Nx+1);
    Matrix M_upd_rho2 = Matrix(Nx+1,Nx+1);
    // For the Newton_Part for v
    Vector T = Vector(Nx+1); Matrix Jac_T = Matrix(Nx+1,Nx+1); Vector v = Vector(Nx+1); Vector delta = Vector(Nx+1); Vector old_v = Vector(Nx+1); double h = 1E-10; int it_newt(0); T[0] = 1.0;
    while( ( u - old_u).linfty()/u.linfty() > 1E-7 && iter < itermax)
    {
        
        // Re-initialize the Matrix
        for(int i= 0 ; i <= Nx ; i++)
        {
            for(int j = 0; j <= Nx ; j++)
            {
                M_upd_rho[i][j] = 0.;
                Jac_T[i][j] = 0;
                
            }
            
        }
        // Update rho with the new computed velocity
       for(int i = 0; i <= Nx ; i++)
        {
            ir = mod(i+1,Nx+1); il = mod(i-1,Nx+1);
            M_upd_rho[i][i]    = 1  + (dt/dx) * (g(u[i])  + g(u[il]) - u[il] );
            M_upd_rho[i][ir]   = -(dt/dx)   * (g(u[i]) - u[i]);
            M_upd_rho[i][il]   = -(dt/dx)   * g(u[il]);
            P[i] = i ; // Initialize P to identity
        }
        ok = Gauss(M_upd_rho,P);
        if(ok != 0){rho = Solve(M_upd_rho,P,rho_at_step_n);} // Solver pas assez robuste ?
        
        
        
        // udpate_phi
        //Compute the new potential
        solve_non_linear_Poisson(rho_tot, phi);
        
        //udpate u using a Newton algorithm
        // Initialize Newton
        v = u ;
        while ( T.linfty() > 1E-10 && it_newt < itermax)
        {
            for(int i = 0 ; i <= Nx; i++)
            {
                irr = mod(i+2,Nx+1); ir = mod(i+1,Nx+1); il = mod(i-1,Nx+1); ill = mod(i-2,Nx+1);
                //Compute Flux of momentum
                Q_r = 0.5 * (G(rho[ir],rho[irr],u[ir]) + G(rho[i],rho[ir],u[i]));   // Q_{i+1}
                Q_l = 0.5 * (G(rho[i],rho[ir],u[i]) + G(rho[il],rho[i],u[il])); // Q_{i}
                
                //Compute the Jacobian of the Non linear Application T  for which we look for v such that T(v)= 0.
                
                Jac_T[i][i]  = 0.5 * (rho[ir]+rho[i]) + (dt/dx) * (max(Q_r,0.)  + max(Q_l,0.) - Q_l) -  dt * ((phi[ir]-phi[i])/dx) * (t_rho(rho[i],rho[ir],v[i]+h) - t_rho(rho[i],rho[ir],v[i]-h))/(2*h);
                Jac_T[i][ir] = -(dt/dx)   * (max(Q_r,0.) - Q_r);
                Jac_T[i][il] = -(dt/dx)   * max(Q_l,0.);
                
                // Compute T(v)
                T[i] = 0.5 * (rho[ir] + rho[i]) * v[i] - 0.5 * (rho_at_step_n[ir] + rho_at_step_n[i]) * u_at_step_n[i] + (dt/dx) * (max(Q_r,0.)  + max(Q_l,0.) - Q_l) * v[i] -(dt/dx)   * (max(Q_r,0.) - Q_r)*v[ir] -(dt/dx) * max(Q_l,0.) * v[il] - dt * ((phi[ir]-phi[i])/dx) * t_rho(rho[i],rho[ir],v[i]);
                // Permutation Matrix for LU Factorization
                P[i] = i;
             
            }
            
            ok = Gauss(Jac_T,P);
            if(ok != 0){delta = Solve(Jac_T,P,T);}
            v = v - delta;
            
            //cout << "                   internal Newton_iteration for the velocity  : " << it_newt << "    T(v(u)^{k}) = " << T.linfty() << endl;
            it_newt++;
        }
        old_u = u;
        u = v ; //old_u = u;
        //cout << "       internal fixed-point iteration = " << iter << " ||u^{"<<iter+1<<"} - u^{" << iter << "} || / ||u^{ "<< iter <<"}||   = " <<(u - old_u).linfty()/old_u.linfty() << endl;
        iter++;
    }
    
    //cout <<"    EndCall : solve_pressure-less_Euler-Poisson." << endl;
}



void plot_sol(int iter, Matrix & rho, Matrix & u, Vector & phi)
{
  //Plot the sol over 1-period
  Vector source_term = Vector(Nx+1);
  string file_name_f = "sim/f_at_"+to_string(iter)+".dat";
  string file_name_E = "sim/E_at_"+to_string(iter)+".dat";
  ofstream file_f(file_name_f); ofstream file_E(file_name_E);
    double x_i; double alpha(0); double omega(0); double f(0);
    int ir; int il;
    for(int i = 0 ; i<=  Nx ; i++)
    {
        ir = mod(i+1,Nx+1); il = mod(i-1,Nx+1);
        x_i = i * dx;
        file_E << x_i << "  " << (phi[ir]-phi[i])/dx << endl;
        f = 0;
        for(int j= 0 ; j < nG ; j++)
        {
            if(hermite_quad == true)
            {
                alpha = quad._X[j];
                omega = quad._W[j]  * mean_f0(alpha) * exp(alpha*alpha);
            }
            else
            {
                alpha = vmin + j * ( (vmax-vmin)/(nG-1) );
                omega = mean_f0(alpha)/SF0;
            }
            
            //source_term[i] = t_rho(rho[i],rho[ir],u[i]) * ((phi[ir]-phi[i])/dx);
            file_f << x_i << "  " << u[j][i] << "  " <<   rho[j][i] <<  endl;
        }
        file_f << endl;
    }
    file_f.close(); file_E.close();
}

