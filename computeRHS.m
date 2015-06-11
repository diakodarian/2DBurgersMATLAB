% function [getRHS] = computeRHS(u,nu,jacobian,Nx,Ny,Lx)
% 
%
%
%   Author:   Diako Darian
%   Date:     11.07.2015
%
%
% This function computes the RHS of 2D Burger's equation
% 
%         du = -uu_x-vu_y + nu (u_xx+u_yy) (+ MMS-term),
%
%
%      where  u = (u,v) is a 2D vector (Velocity field)
%      
%      BC for u: periodic in x-direction
%                and Neumann BC in y-direction u_y(x,-1) = u_y(x,1) = 0   
% 
%      BC for v: periodic in x-direction
%                and Dirichlet BC in y-direction v(x,-1) = v(x,1) = 0 
% 
%   Purpose    :  RHS of the 2D Burgers' equation
%    
%                        u_t = -(u.grad)u + nu (del^2)u
%
%                 is solved by Chebyshev spectral collocation method 
%                 using FFT in y-direction, and by pseudo-spectral
%                 Fourier-Galerkin method in x-direction.
%
% 
%
%----------------------------ooooooooo-------------------------------------

function [getRHS] = computeRHS(u,nu,jacobian,Nx,Ny,Lx)


%--------------------------------------------------------------------------
%                         Nonlinear term:
%--------------------------------------------------------------------------

% Filter for dealiasing nonlineat convection:
k = [0:Nx/2-1 0 -Nx/2+1:-1];
kmax = (2/3)*(Nx/2 + 1);
dealias = abs(k)>kmax;

uu_x = zeros(Ny+1,Nx); uv_x = zeros(Ny+1,Nx);
vu_y = zeros(Ny+1,Nx); vv_y = zeros(Ny+1,Nx);

Du1 = zeros((Ny+1)*Nx,2);

uux = reshape(u(:,1), [Ny+1, Nx]);
uuy = reshape(u(:,2), [Ny+1, Nx]);

% Fourier transform in x-direction (Periodic direction)
% u d_x u
for i = 1:Ny+1                % 1nd derivs wrt x in each row
  vv = uux(i,:);
  v_hat = fft(vv); vvr = ifft(v_hat); vvr = vvr.^2;
  vvr_hat = fft(vvr);
  w_hat = (pi/Lx)*1i*[0:Nx/2-1 0 -Nx/2+1:-1].* vvr_hat;
  w_hat(dealias) = 0;
  w = real(ifft(w_hat)); 
  uu_x(i,:)= w;
end
% u d_x v
for i = 1:Ny+1                % 1nd derivs wrt x in each row
  vvy = uux(i,:); vvx = uux(i,:);
  vy_hat = fft(vvy); vvyr = ifft(vy_hat);
  vx_hat = fft(vvx); vvxr = ifft(vx_hat);
  vvr = vvxr.* vvyr;
  vvr_hat = fft(vvr);
  w_hat = (pi/Lx)*1i*[0:Nx/2-1 0 -Nx/2+1:-1].* vvr_hat;
  w_hat(dealias) = 0;
  w = real(ifft(w_hat)); 
  uv_x(i,:)= w;
end

% Chebyshev in y-direction (non-periodic)
% v d_y v
for j = 1:Nx                % 1nd derivs wrt y in each column
    uy = uuy(:,j); 
    vv_y(:,j) = jacobian*chebfft(uy)'; 

    N = length(uy)-1; 
    v_hat = fft(vv_y(:,j)); uy_hat = fft(uy);  
   
    vr = ifft(v_hat); uyr = ifft(uy_hat);
    
    v2 = vr.*uyr; v2_hat = fft(v2);
    % Filter for dealiasing nonlineat convection:
    k = [0:N/2-1 0 -N/2+1:-1];
    kmax = (2/3)*(N/2 + 1);
    dealias = abs(k)>kmax;
    v2_hat(dealias) = 0;
    vdv = real(ifft(v2_hat));
    vv_y(:,j) = vdv;
end 

% v d_y u
for j = 1:Nx                % 1nd derivs wrt y in each column
    uy = uuy(:,j); ux = uux(:,j); 
    vu_y(:,j) = jacobian*chebfft(ux)'; 
    % Impose Neumann BC:
    vu_y(1,j) = 0; vu_y(end,j) = 0;

    N = length(uy)-1; 
    v_hat = fft(vu_y(:,j)); uy_hat = fft(uy);  
   
    vr = ifft(v_hat); uyr = ifft(uy_hat);
    
    v2 = vr.*uyr; v2_hat = fft(v2);
    % Filter for dealiasing nonlineat convection:
    k = [0:N/2-1 0 -N/2+1:-1];
    kmax = (2/3)*(N/2 + 1);
    dealias = abs(k)>kmax;
    v2_hat(dealias) = 0;
    vdv = real(ifft(v2_hat));
    vu_y(:,j) = vdv;
end 

du_x = -uu_x - vu_y; du_y = -uv_x - vv_y;
Du1(:,1) = du_x(:); Du1(:,2) = du_y(:);
%-------------------------------------------------------------------------



%--------------------------------------------------------------------------
%                         Diffusion term:
%--------------------------------------------------------------------------
Du = zeros((Ny+1)*Nx,2);

for k = 1:2
    uxx = zeros(Ny+1,Nx); uyy = zeros(Ny+1,Nx);
    uu = reshape(u(:,k), [Ny+1, Nx]);
    
    % Fourier transform in x-direction (Periodic direction)
    for i = 1:Ny+1                % 2nd derivs wrt x in each row
      vv = uu(i,:);
      v_hat = fft(vv);
      w_hat = (pi/Lx)^2*((1i*[0:Nx/2-1 0 -Nx/2+1:-1]).^2).* v_hat;
      w = real(ifft(w_hat)); 
      uxx(i,:)= w;
    end

    % Chebyshev in y-direction (non-periodic)
    for j = 1:Nx                % 2nd derivs wrt y in each column
      uy = uu(:,j); 
      uyy(:,j) = jacobian*chebfft(uy)'; 
      % Impose Neumann BC:
      if k == 1
          uyy(1,j) = 0; uyy(end,j) = 0;
      end
      uyy(:,j) = jacobian*chebfft(uyy(:,j))';
    end 
    du = uxx + uyy;
    Du(:,k) = du(:);
end

%--------------------------------------------------------------------------
%                       Assemble RHS:
%--------------------------------------------------------------------------

getRHS = Du1 + nu*Du; 