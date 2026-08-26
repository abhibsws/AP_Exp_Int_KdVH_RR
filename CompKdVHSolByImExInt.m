function [ImExTime,dt,u,v,w] = CompKdVHSolByImExInt(Len, tau, ImExInt, type, tf, N, u0, v0, w0, dt)
%=========================================================================%
%  CompKdVHSolByImExInt: solves the KdVH system using ImEx methods        %
%-------------------------------------------------------------------------%
%  Inputs:                                                                %
%    Len      - length of the computational domain                        %
%    tau      - relaxation parameter                                      %
%    ImExInt  - name of implicit-explicit integrator                      %
%    type     - 'I' or 'II'                                               %
%    tf       - Final simulation time                                     %
%    N        - Number of spatial grid points                             %
%    u0, v0, w0 - Initial conditions for solution components              %
%    dt       - time step                                                 %
%                                                                         %
%  Outputs:                                                               %
%    ImExTime  - runtime                                                  %
%    u   - Solution component u at final time                             %
%    v   - Solution component v at final time                             %
%    w   - Solution component w at final time                             %
%=========================================================================%
    % Butcher tableau for the ImEx method
    [rkim,bim,cim, rkex, bex, cex] = ImExButcherTableau(type,ImExInt);    
    % Wave numbers
    xi = (-N/2:N/2-1);        % Indices of the wave numbers
    xi = fftshift(xi);      % Reordering of the indices
    xi = xi*2*pi/Len;         % Scaling by 2π/L. 
    
    % Filtering
    filtr = ones(size(xi));  % Create an array of ones with the same size as xi
    xi_max = max(abs(xi));   % Compute the maximum absolute value of xi


    filtr(abs(xi) > xi_max * (2/3)) = 0;  % Set elements to 0 where condition is met
    
    % The initial profile, obtained from petviashvili reference solution
    u = u0; v = v0; w = w0;
    % Time step 
    nt = round(tf/dt);
    dt = tf/nt; tn = 0;
 
    Y = [u;v;w]; % u,v,w are row vectors, each of size 1 by N, so Y is 3 by N.
     
    
    % construct the block diagonal matrix for the lienar part of RHS
    L_tilde = cell(1, N); 
    for i = 1:N
        L_tilde{1,i} = [  0,  0,  -1j*xi(i);
                          0,  1j*xi(i)/tau, -1/tau;
                         -1j*xi(i)/tau,  1/tau,  0];
    end 
    L = sparse(blkdiag(L_tilde{:}));
    
    s = size(rkex,1);
    
    % Precompute LU factorization for (I - dt * rkim(j,j) * L)
    LU_factors = cell(1, s);
    I_sparse = speye(3*N);  % Identity matrix of size 3N x 3N
    for j = 1:s
        A_j = I_sparse - dt * rkim(j,j) * L; % System matrix
        [L_j, U_j, P_j] = lu(A_j); % LU factorization
        LU_factors{j} = {L_j, U_j, P_j}; % Store factors
    end 
    %
    tic          
    % time loop: Solving the ODE system in spatial domain
    for i = 1:nt
        g = zeros(3*N,s);     % for storing intermediate stages as column vectors
        Rex = zeros(3*N,s);       % Evauation of nonlinear part in the rhs at g_j
        Rim = zeros(3*N,s);          % Evauation of linear part in the rhs at g_j
        Yh = reshape(fft(Y,[],2),3*N,1);
        for j = 1:s
            rhs = Yh;
            if j>1
                rhs = rhs + dt*Rex(:,1:j-1)*(rkex(j,1:j-1))' + dt*Rim(:,1:j-1)*(rkim(j,1:j-1))';
            end 
            % linear solve for the stages: LU decomposition appraoch
            L_j = LU_factors{j}{1}; U_j = LU_factors{j}{2}; P_j = LU_factors{j}{3};
            ysol = L_j \ (P_j * rhs);  % Forward substitution
            gh_j = U_j \ ysol;          % Backward substitution
    
            % linear solve for the stages: Maltab backslash
            %gh_j = (I_sparse - dt * rkim(j,j) * L) \ rhs;
    
            % Calculating the right hand side
            Rex(:,j) = nonlin_h_fun(gh_j,xi,filtr);
            Rim(:,j) = L*gh_j;
            g(:,j) = gh_j;
        end
        Yh_new = Yh + dt*(Rex*bex+Rim*bim); % solution update in the Fourier space
        Y_new = real(ifft(reshape(Yh_new, 3, N),[],2)); % This applies inverse fft row-wise
        Y = Y_new; tn = tn+dt;    
    end
    u = Y(1,:); v = Y(2,:); w = Y(3,:);
    ImExTime = toc;
end
%------------------------------functions----------------------------------%
% Nonlinear part of the RHS
function Yh_rhs = nonlin_h_fun(Yh,xi,filtr)
    Yh_rhs = zeros(length(Yh),1);  % Initialize Y with zeros
    uhat = Yh(1:3:end,1); 
    uh_rhs = -(1j*xi'/2).*fft(real(ifft(uhat.*filtr')).^2);
    Yh_rhs(1:3:end) = uh_rhs;
end
%--------------------------------End--------------------------------------%
