
% function [tdata,udata] = PseudoSpectralFourier(Nx,Ny,Lx,...,u0,tmax)
% 
%
%  This function advances RHS of the 2D Burger's equation in time
%  by means of fourth order Runge-Kutta method, 
%  in the case of pseudo-spectral Fourier-Galerkin method.
%
% 
%      u_t = -uu_x-vu_y + nu (u_xx+u_yy) + MMS-term
%
%         
%
% Author:   Diako Darian
% Date:     11.07.2015
% 
% 
% 
% Purpose    :  Time integration of 2D Burgers' equation
%    
%                        u_t = -(u.grad)u + nu (del^2)u
%
% by RK4.
%
%----------------------------ooooooooo-------------------------------------

function [tdata,udata] = PseudoSpectralFourier(Nx,Ny,Lx,Ly,nu,x0,y0,u0,tmax)

% Set up grid and initial data:
  dt = 0.1/(Nx^2+Ny^2); hx = 2*pi/Nx; hy = 2*pi/Ny;
  x = x0; y = y0;
  
  h = dt;                                        % step size
  % RK4 parameters:
  a = [1/6 1/3 1/3 1/6]*h;
  b = [.5 .5 1]*h;
  
% Initial condition:
  clf, drawnow, set(gcf,'renderer','zbuffer')
  u = u0;
  v = fft2(u); 
  

% Solve PDE and plot results:
  clf, plotgap = round(tmax/dt);
  [xx,yy] = meshgrid(x,y);

  for n = 0:2*plotgap
    t = n*dt; 
    % Plot the results:
    if rem(n+.5,plotgap)<1
        vv = real(ifft2(v)); 
        vvv = reshape(vv(:,1), [Ny, Nx]);
        subplot(3,1,n/plotgap+1), mesh(xx,yy,vvv), view(-10,60)
        axis([-Lx Lx -Ly Ly -0.15 1]), colormap([.4 .3 .2]);
        text(-2.5,1,.5,['t = ' num2str(t)],'fontsize',18),
        set(gca,'ztick',[]),xlabel x, ylabel y, zlabel u, grid off, drawnow
    end
    % RK4: 
    v1 = v; v0 = v;
    for rk = 1:4
        du = computeRHS_Pseudo(v,nu,Nx,Ny,Lx,Ly,xx,yy,t);
        if rk < 4
            v = v0; v = v + b(rk)*du;
        end
        v1 = v1 + a(rk)*du;
    end
    v = v1; 
  end

%------------------------ooooooooooooo-------------------------------------
%                           MMS Test
%--------------------------------------------------------------------------
figure(2);
for n = 0:2*plotgap
    t = n*dt; 

    ux = sin(xx).*cos(yy)*exp(-2*nu*t);
    uy = -cos(xx).*sin(yy)*exp(-2*nu*t);
    
    if rem(n+.5,plotgap)<1
        subplot(3,1,n/plotgap+1), mesh(xx,yy,ux), view(-10,60)
        axis([-Lx Lx -Ly Ly -0.15 1]), colormap([.4 .3 .2]);
        text(-2.5,1,.5,['t = ' num2str(t)],'fontsize',18),
        set(gca,'ztick',[]),xlabel x, ylabel y, zlabel u, grid off, drawnow
    end
end
%------------------------ooooooooooooo-------------------------------------