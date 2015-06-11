%  function [getRHS] = computeRHS_Pseudo(v, nu, Nx, Ny, Lx, Ly,x,y,t)
% 
%
%  This function computes the RHS of 2D Burger's equation in the case of
%  pseudo-spectral Fourier-Galerkin method.
%
% 
%  Compute RHS: du = -uu_x-vu_y + nu (u_xx+u_yy) + MMS-term
%
%         
%
% Author:   Diako Darian
% Date:     11.07.2015
% 
% 
% 
% Purpose    : computeRHS.m computes RHS of 2D Burgers' equation
%    
%                        u_t + (u.grad)u = nu (del^2)u
%
% by FFT.
%
%----------------------------ooooooooo-------------------------------------

function [getRHS] = computeRHS_Pseudo(v, nu, Nx, Ny, Lx, Ly,x,y,t)


kx = [0:Nx/2-1 0 -Nx/2+1:-1]*(pi/Lx); 
ky = [0:Ny/2-1 0 -Ny/2+1:-1]*(pi/Ly); 
  
[Kx,Ky] = meshgrid(kx,ky); 
K = [Kx(:), Ky(:)]; K2 = K.^2;

%--------------------------------------------------------------------------
%        Filter for dealiasing nonlineat convection:
%--------------------------------------------------------------------------

kxmax = (2/3)*(Nx/2 + 1)*(pi/Lx);
kymax = (2/3)*(Ny/2 + 1)*(pi/Ly);
dealias = zeros(size(K));
dealias(:,1) = abs(K(:,1))>kxmax;
dealias(:,2) = abs(K(:,2))>kymax;

%--------------------------------------------------------------------------
%                         Nonlinear term
%--------------------------------------------------------------------------

u = ifft2(v);
u2 = u.*u;

u2_hat = fft2(u2);
u_nonlinear = -1i*K.*u2_hat;

% Filter for dealiasing nonlineat convection:

u_nonlinear = u_nonlinear.*dealias;


%--------------------------------------------------------------------------
%                          MMS Test
%--------------------------------------------------------------------------

fx = -sin(x).*cos(x)*exp(-4*nu*t);
fy = -sin(y).*cos(y)*exp(-4*nu*t);

fx = fx(:); fy = fy(:);

F = zeros(size(K));

F(:,1) = fx; F(:,2) = fy; 

F_hat = fft2(F);
%----------------------oooooooooooooo--------------------------------------

%--------------------------------------------------------------------------
%                       Assemble the RHS:
%--------------------------------------------------------------------------

getRHS = u_nonlinear - nu*K2.*v + F_hat;  
