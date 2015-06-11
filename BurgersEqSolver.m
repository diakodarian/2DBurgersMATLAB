% 
% BurgersEqSolver.m 
%         
%
% Author:   Diako Darian
% Date:     11.07.2015
% 
% 
% 
% Purpose    : BurgersEqSolver solves 2D Burgers' equation
%    
%                       u_t = -(u.grad)u + nu (del^2)u
%
% This equation will be solved by using Chebyshev spectral collocation
% method in the case of non-periodic BC, and by a pseudo-spectral 
% Fourier-Galerkin method in the case of a periodic BC.
%
%
% This code solves the above equation for mixed boundary conditions 
% Choose BC = 0 for 2D periodic BC and Fourier-Galerkin method
% Choose BC = 1 for one periodic and one non-periodic BC
%
%     --------------oooooo---------------------
%
% MMS test problem for Fourier-Galerkin method: 
%   
%           Initial condition   u = sin(x).*cos(y);
%                               v = -cos(x).*sin(y);
%           
%           Periodic BC 
% 
%          
%
%     --------------oooooo---------------------
%
% Test problem for Chebyshev spectral collocation method: 
%   
%           Initial condition u(x,0) = u0 = 2x 
%           
%           Dirichlet BC u(-1,t) = -2/(1+2t) and u(1,t) = 2/(1+2t)
% 
%           This problem has analytical solution u(x,t) = 2x/(1+2t) 

BC = 1;

%--------------------------------------------------------------------------
%                   Pseudo-Spectral Fourier-Galerkin Method
%--------------------------------------------------------------------------
if BC == 0

 % Grid and initial conditions:
  Nx = 32; Ny = 32; nu = 0.1; Lx = pi; Ly = pi; 
  hx = 2*pi/Nx; x = hx*(1:Nx); x = Lx*(x-pi)/pi;
  hy = 2*pi/Ny; y = hy*(1:Ny); y = Ly*(y-pi)/pi;
  [X,Y] = meshgrid(x,y);
  ux = sin(X).*cos(Y);
  uy = -cos(X).*sin(Y);
  u0 = [ux(:), uy(:)]; 
  tmax = 0.1; 
  PseudoSpectralFourier(Nx,Ny,Lx, Ly,nu,x,y,u0,tmax)
 
%--------------------------------------------------------------------------
%                  Chebyshev Spectral Collocation Method
%--------------------------------------------------------------------------
elseif BC == 1

  % Grid and initial conditions:
  Nx = 36; Ny = 20; nu = 0.1; Lx = 3;
  a = -10; b = 10; Ly = b - a; 
  hx = 2*Lx/Nx; x = -Lx + hx*(1:Nx);
  y = cos(pi*(0:Ny)/Ny); Jacobian = 2/Ly;
   
  [X,Y] = meshgrid(x,y);
  ux = exp(-8*((X+1.5).^2+Y.^2));
  uy = exp(-8*((X+1.5).^2+Y.^2));
  %uu = [ux, uy];
  u0 = [ux(:), uy(:)];
  tmax = 0.2;
  SpectralChebyshevFFT(Nx,Ny,a,b,Lx, Ly,nu,x,y,u0,tmax,Jacobian)
  
end