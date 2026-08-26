%--------------------------------------------------------------------------
% This script computes and saves reference solutions for the KdVH system 
% using the Petviashvili iteration method. These solutions will be used for 
% further numerical studies.
%--------------------------------------------------------------------------
%compute reference solution using Petviashvili
xL = -40; xR = 40; tau =1e-8; c = 1.2; Nref = 2^10; maxiter = 1000; tol = 5.0e-12;
[uP, vP, wP, xP] = petviashvili_kdvh(xL,xR,tau,c,Nref,maxiter,tol);

filename = sprintf('RefSolPetviashvili_c_%.1f_tau_%.1e.mat', c, tau); 
save(filename, 'uP', 'vP', 'wP', 'xP');

% To load the file whereever required
% filename = sprintf('RefSolPetviashvili_c_%.1f_tau_%.e.mat', c, tau);
% load(filename);
%-------------------------------Functions---------------------------------%
function [u, v, w, x] = petviashvili_kdvh(xmin,xmax,tau,c,N,maxiter,tol)
%-------------------------------------------------------------------------%
% Function: petviashvili_kdvh                                             
% This function computes a stationary solution of the KdVH system using the 
% Petviashvili iteration method.                                                                                 
% Inputs:                                                                 
% - xmin   : Left boundary of the spatial domain                         
% - xmax   : Right boundary of the spatial domain                        
% - tau    : Small parameter in the KdVH equation                        
% - c      : Wave speed parameter                                        
% - N      : Number of spatial grid points                               
% - maxiter: Maximum number of iterations for convergence                
% - tol    : Tolerance for convergence criterion                         
%                                                                         
% Outputs:                                                                
% - u : Computed stationary solution                                     
% - v : First derivative of u                                     
% - w : Second derivative of u                                         
% - x : Spatial grid points                                              
%-------------------------------------------------------------------------%
    % Construct the Fourier differentiation matrix
    [D,ix] = fourier_first_derivative_matrix(xmin, xmax, N); x = ix;
    % Define the linear operator L 
    L = c/((1+c*tau)*(1-c^2*tau))*eye(N) - D^2;
    % Define the nonlinear function N(u) 
    N_fun = @(u) 0.5/((1+c*tau)*(1-c^2*tau)) * u.^2 + (c*tau)/(1-c^2*tau)*D*(u.*(D*u));
    % Parameter gamma for Petviashvili method
    gamma = 2;
    % Spatial grid
   
    % Initial condition
    u = exp(-ix.^2);
    % Petviashvili iteration
    for iter = 1:maxiter
        res = norm(L*u-N_fun(u), inf);
        if res < tol
            fprintf('Converged in %d iterations with residual %e.\n', iter, res);
            break;
        end
        % Compute the stabilizing factor m(u)
        Lu = L * u;
        Nu = N_fun(u);
        m = (sum(Lu.*u)/sum(Nu.*u));
        %m = integrate(u .* Lu, D) / integrate(u .* Nu, D);
        % Update the solution
        u = m^gamma*(L\N_fun(u));
    end
    w = c*u-u.^2/2; v = D*u - c*tau*D*w;
    disp(['Final residual: ', num2str(norm(L*u -  N_fun(u))), '.']);
end

%------------------------------functions----------------------------------%
% Fourier differetiation matrix
function [D,ix] = fourier_first_derivative_matrix(xmin, xmax, Nx)
    L = xmax - xmin; 
    x = linspace(xmin, xmax, Nx+1)'; ix = x(1:end-1); Nf = Nx;
    
    % To compute fhat
    for ii = 1 : Nf
        kk = ii - Nf / 2;
        for jj = 1 : Nx
            omega(ii, jj) = exp(-1i*kk*ix(jj) * 2 * pi / L) / Nx;
        end
    end
    
    % Now take the inverse using fhat
    for ii = 1 : Nx
        for jj = 1 : Nf
            kk = jj - Nf / 2;
            iomega(ii, jj) = 1i * kk  * 2 * pi / L * exp(1i * kk * ix(ii) * 2 * pi / L);
        end
    end
    
    % Cheating ;-)
    D = real(iomega * omega);
end
%--------------------------------End--------------------------------------%