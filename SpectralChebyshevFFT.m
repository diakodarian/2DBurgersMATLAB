% function [tdata,udata] =
% SpectralChebyshevFFT(Nx,Ny,a,b,nu,x0,y0,u0,tmax,jacobian)
% 
%
%
%   Author:   Diako Darian
%   Date:     11.07.2015
%
%
%  This function advances RHS of the 2D Burger's equation in time
%  by means of fourth order Runge-Kutta method, 
%  in the case of Chebyshev spectral collocation method using FFT.
%
% 
%      u_t = -uu_x-vu_y + nu (u_xx+u_yy) + MMS-term,
%
%      where  u = (u,v) is a 2D vector (Velocity field)
%      
%      BC for u: periodic in x-direction
%                and Neumann BC in y-direction u_y(x,-1) = u_y(x,1) = 0   
% 
%      BC for v: periodic in x-direction
%                and Dirichlet BC in y-direction v(x,-1) = v(x,1) = 0 
% 
%   Purpose    :  Time integration of 2D Burgers' equation
%    
%                        u_t = -(u.grad)u + nu (del^2)u
%
%                 by RK4.
%
%----------------------------ooooooooo-------------------------------------

function [tdata,udata] = SpectralChebyshevFFT(Nx,Ny,a,b,Lx,Ly,nu,x0,y0,u0,tmax,jacobian)

% Set up grid and initial data:
    dt = 1/(Nx+Ny^2); x = x0; y = y0; 
    yi = .5*b*(1+y) + .5*a*(1-y);

% Initial condition:
    clf, drawnow, set(gcf,'renderer','zbuffer')
    u = u0;

% Solve 2D Burgers' equation and plot results:
    clf, plotgap = round(tmax/dt);
    [xx,yy] = meshgrid(x,y); BCa = yy(:) == a; BCb = yy(:) == b; 
        
%--------------------------------------------------------------------------
%           Time-stepping using RK4 method
%--------------------------------------------------------------------------

    % RK4 parameters:
    a_rk4 = [1/6 1/3 1/3 1/6]*dt;
    b_rk4 = [.5 .5 1]*dt;

    for n = 0:2*plotgap
        t = n*dt; 

        % Plot the results:
        if rem(n+.5,plotgap)<1
            vv = reshape(u(:,1), [Ny+1, Nx]);
            subplot(3,1,n/plotgap+1), mesh(xx,yy,vv), view(-10,60)
            %axis([-Lx Lx a b -0.15 1]), colormap([.4 .3 .2]);
            text(-2.5,1,.5,['t = ' num2str(t)],'fontsize',18),
            set(gca,'ztick',[]),
            xlabel x, ylabel y, zlabel u, grid off, drawnow
        end % end if

        % RK4 loop 
        u1 = u; u0 = u;
        for rk = 1:4
            du = computeRHS(u,nu,jacobian,Nx,Ny,Lx);
            if rk < 4
                u = u0; u = u + b_rk4(rk)*du;
            end
            u1 = u1 + a_rk4(rk)*du;
        end % end for RK4 loop
        u = u1;
        % Impose Dirichlet BC on v in y-direction:
        u(BCa,2) = 0; u(BCb,2) = 0;

    end % end for (time integration)

 %             -------------oooooooo----------------   
                                         

