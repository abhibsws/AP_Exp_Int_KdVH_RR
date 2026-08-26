function [ExpIntTime,dt,u,v,w] = CompKdVHSolByExpInt(Len, tau, exp_int, tf, N, u0, v0, w0, dt)
%=========================================================================%
% CompKdVHSolByExpInt: Solves the KdVH system using exponential methods   %
%-------------------------------------------------------------------------%
%  Inputs:                                                                %
%    tau      - relaxation parameter                                      %
%    CFL      - CFL number for time step calculation                      %
%    exp_int  - name of exponential integrator                            %
%    tf       - Final simulation time                                     %
%    N        - Number of spatial grid points                             %
%    u0, v0, w0 - Initial conditions for solution components              %
%                                                                         %
%  Outputs:                                                               %
%    dt  - Computed time step                                             %
%    x   - Spatial grid points                                            %
%    u   - Solution component u at final time                             %
%    v   - Solution component v at final time                             %
%    w   - Solution component w at final time                             %
%=========================================================================%
    % Wave numbers
    xi = -N/2:N/2-1;        % Indices of the wave numbers
    xi = fftshift(xi);      % Reordering of the indices
    xi = xi*2*pi/Len;         % Scaling by 2π/L
    
    % Filtering
    filtr = ones(size(xi));  % Create an array of ones with the same size as xi
    xi_max = max(abs(xi));   % Compute the maximum absolute value of xi
    filtr(abs(xi) > xi_max * (3/3)) = 0;
    
    % The initial profile, obtained from petviashvili reference solution
    u = u0; v = v0; w = w0;
    nt = round(tf/dt);
    dt = tf/nt; tn = 0;
    
    % Time stepping 
    Y = [u;v;w]; % u,v,w are row vectors, each of size 1 by N, so Y is 3 by N. 
                % Each column corresponds to a wave number.
    
    % Computing the matrix exponential using Lagrange Sylvester formula 
    phi0 = cell(1, N);  % Preallocate cell array to store blocks
    phi02 = cell(1, N);
    phi1 = cell(1, N);  
    phi12 = cell(1, N);
    phi2 = cell(1, N);
    phi22 = cell(1, N);
    phi3 = cell(1, N);
    phi32 = cell(1, N);
    % For ETD5RKF
    phi2hat = cell(1,N); 
    u11 = cell(1,N); u21 = cell(1,N); u31 = cell(1,N);
    u41 = cell(1,N); u51 = cell(1,N); u61 = cell(1,N);
    for i = 1:N
        L = [  0,  0,  -1j*xi(i);
                  0,  1j*xi(i)/tau, -1/tau;
                 -1j*xi(i)/tau,  1/tau,  0];
    
        [phi_02, phi_12, phi_22, phi_32] = phi_funs_on_mat(dt,L/2);
        % Compute phi02, phi_12, phi_12.
        phi02{i} =  phi_02; phi12{i} =  phi_12; phi22{i} =  phi_22; phi32{i} = phi_32;
    
        % Computing exp using Lagrange Sylvester equation
        [phi_0, phi_1, phi_2, phi_3] = phi_funs_on_mat(dt,L);
        % Compute phi0
        phi0{i} = phi_0;
        % Compute phi1
        phi1{i} =  phi_1;
        % Compute phi2
        phi2{i} =  phi_2;
        % Compute Phi3
        phi3{i} =  phi_3;

        % Compute Phi2_hat = Phi2_3/5
        [~, ~, phi_2_hat,~] = phi_funs_on_mat(dt,3*L/5);
        phi2hat{i} = phi_2_hat;

        %
        c = [0, 2/9, 1/3, 3/4, 1, 5/6];
        [u21{i},~,~,~] = phi_funs_on_mat(dt,c(2)*L); 
        [u31{i},~,~,~] = phi_funs_on_mat(dt,c(3)*L); [u41{i},~,~,~] = phi_funs_on_mat(dt,c(4)*L);
        [u51{i},~,~,~] = phi_funs_on_mat(dt,c(5)*L); [u61{i},~,~,~] = phi_funs_on_mat(dt,c(6)*L);
    
        % Post correction 
        %phi1{1} = zeros(3); phi1{1}(1,1) = 1;
        phi1{1} = eye(3) + dt*L/2;
        % Post correction 
        %phi2{1} = zeros(3); phi2{1}(1,1) = 1/2;
        phi2{1} = eye(3)*(1/2) + dt*L/6;
        % Post correction 
        %phi3{1} = zeros(3); phi3{1}(1,1) = (1/6);
        phi3{1} = (1/6)*eye(3)+dt*L/24;
        % Post correction 
        %phi12{1} = zeros(3); phi12{1}(1,1) = 1;
        phi12{1} = eye(3) + dt*L/4;
        % Post correction 
        phi22{1} = eye(3)*(1/2) + dt*L/12;
        % Post correction for the zero frequency (xi=0)
        phi32{1} = (1/6)*eye(3) + dt*L/48; % Approx for phi3(dt*L/2)
        %
        u11{1} = eye(3); u21{1} = eye(3); u31{1} = eye(3);
        u41{1} = eye(3); u51{1} = eye(3); u61{1} = eye(3);
        phi2hat{1} = eye(3)*(1/2) + dt*L/10; 
    end
    
    % construction of the block diagonal matrix for all the phi functions
    Phi0 = sparse(blkdiag(phi0{:}));
    Phi02 = sparse(blkdiag(phi02{:}));
    Phi1 = sparse(blkdiag(phi1{:}));
    Phi12 = sparse(blkdiag(phi12{:}));
    Phi2 = sparse(blkdiag(phi2{:}));
    Phi22 = sparse(blkdiag(phi22{:}));
    Phi3 = sparse(blkdiag(phi3{:}));
    Phi32 = sparse(blkdiag(phi32{:}));
    % 
    U11 = sparse(blkdiag(u11{:})); U21 = sparse(blkdiag(u21{:}));
    U31 = sparse(blkdiag(u31{:})); U41 = sparse(blkdiag(u41{:}));
    U51 = sparse(blkdiag(u51{:})); U61 = sparse(blkdiag(u61{:}));
    Phi2hat = sparse(blkdiag(phi2hat{:}));
    I_3N = speye(3*N);

    % Precompute ETD5RKF update matrices outside the loop to save time
    if strcmp(exp_int, 'ETD5RKF')
        b11 = (47/150)*Phi1 - (188/75)*Phi2 + (94/15)*Phi3;
        % b12 = zeros
        b13 = (-43/25)*Phi1 + (132/5)*Phi2 - 66*Phi3;
        b14 = (4124/75)*Phi1 - (6152/15)*Phi2 + (2704/3)*Phi3;
        b15 = (189/10)*Phi1 - (662/5)*Phi2 + 284*Phi3;
        b16 = (-1787/25)*Phi1 + (12966/25)*Phi2 - (5628/5)*Phi3;

        % 
        v11 = Phi0;
        
        % Stage coefficients for A(z)
        a21 = -2/3*Phi2 + (10/9)*Phi2hat;
        a31 = (569/11544)*Phi2 + (1355/11544)*Phi2hat;
        a32 = (-831/3848)*Phi2 + (2755/3848)*Phi2hat;
       
        a41 = -(77157/61568)*Phi2 + (143535/61568)*Phi2hat;
        a42 = (587979/61568)*Phi2 - (821745/61568)*Phi2hat;
        a43 = (-405/64)*Phi2 + (675/64)*Phi2hat;
        
        a51 = (655263/7696)*Phi2 - (2031205/23088)*Phi2hat;
        a52 = (-1148769/7696)*Phi2 + (1252665/7696)*Phi2hat;
        a53 = (1593/40)*Phi2 - (405/8)*Phi2hat;
        a54 = (144/5)*Phi2 - (80/3)*Phi2hat;

        a61 = -(2212835/277056)*Phi2 + (6888625/831168)*Phi2hat;
        a62 = (477285/30784)*Phi2 - (496525/30784)*Phi2hat;
        a63 = (-39/16)*Phi2 + (65/16)*Phi2hat;
        a64 = (-4/9)*Phi2 + (20/27)*Phi2hat;
        a65 = (-185/96)*Phi2 + (575/288)*Phi2hat;
    end
    %
    tic
    % Time stepping
    for n = 1:nt
        Yh = reshape(fft(Y,[],2),3*N,1);
        % Different exponential tme integrators
        switch exp_int
            %% Lawson methods
            case 'Lawson-Euler' % Lawson-Euler method of order 1 and stiff order 1
                Yh_rhs = nonlin_h_fun(Yh,xi,filtr);
                Yh_new = Phi0*Yh + dt*Phi0*Yh_rhs;
            case 'Lawson2b' % Lawson method of order 2 and stiff order 1
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi0*Yh + dt*Phi0*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (1/2)*Phi0*Yh1_rhs + (1/2)*Yh2_rhs );         
            case 'Lawson4' % Lawson method of order 4 and stiff order 1
                % 1st stage
                %Yh1 = Yh;
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi02*Yh + (dt/2)*Phi02*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                Yh3 = Phi02*Yh + (dt/2)*Yh2_rhs;
                Yh3_rhs = nonlin_h_fun(Yh3,xi,filtr); % function eval in frequency space
                % 4th stage
                Yh4 = Phi0*Yh + dt*Phi02*Yh3_rhs;
                Yh4_rhs = nonlin_h_fun(Yh4,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (1/6)*Phi0*Yh1_rhs + (1/3)*Phi02*Yh2_rhs +...
                    (1/3)*Phi02*Yh3_rhs + (1/6)*Yh4_rhs );
            %% ETD (Exponential Time Differencing) methods
            case 'Norsett-Euler' % Norsett-Euler method of order 1 and stiff order 1
                Yh_rhs = nonlin_h_fun(Yh,xi,filtr);
                Yh_new = Phi0*Yh + dt*Phi1*Yh_rhs;
            case 'ETD2RK' % ETD2RK method of order 2 and stiff order 2
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi0*Yh + dt*Phi1*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (Phi1 - Phi2)*Yh1_rhs + Phi2*Yh2_rhs );
            case 'ETD3RK' % ETD3RK method of order 3 and stiff order 2
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi02*Yh + dt*(1/2)*Phi12*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                Yh3 = Phi0*Yh + dt*( -Phi1*Yh1_rhs + 2*Phi1*Yh2_rhs );
                Yh3_rhs = nonlin_h_fun(Yh3,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (Phi1-3*Phi2+4*Phi3)*Yh1_rhs + ...
                                (4*Phi2-8*Phi3)*Yh2_rhs +...
                                (-Phi2+4*Phi3)*Yh3_rhs );
            case 'ETD4RK' % ETD4RK method of order 4 and stiff order 2
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi02*Yh + dt*(1/2)*Phi12*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                Yh3 = Phi02*Yh + dt*(1/2)*Phi12*Yh2_rhs;
                Yh3_rhs = nonlin_h_fun(Yh3,xi,filtr); % function eval in frequency space
                % 4th stage
                Yh4 = Phi0*Yh + dt*( (1/2)*Phi12*(Phi02-I_3N)*Yh1_rhs + Phi12*Yh3_rhs );
                Yh4_rhs = nonlin_h_fun(Yh4,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (Phi1-3*Phi2+4*Phi3)*Yh1_rhs + ...
                                (2*Phi2-4*Phi3)*Yh2_rhs +...
                                (2*Phi2-4*Phi3)*Yh3_rhs +...
                                (-Phi2+4*Phi3)*Yh4_rhs );
            case 'ETD5RKF' % 6-stage, 5th order with stiff order 1
                % 1st stage
                Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh1, xi, filtr);            
                % 2nd stage
                Yh2 = U21*Yh + dt*( a21*Yh1_rhs );
                Yh2_rhs = nonlin_h_fun(Yh2, xi, filtr);                
                % 3rd stage
                Yh3 = U31*Yh + dt*( a31*Yh1_rhs + a32*Yh2_rhs );
                Yh3_rhs = nonlin_h_fun(Yh3, xi, filtr);                
                % 4th stage
                Yh4 = U41*Yh + dt*( a41*Yh1_rhs + a42*Yh2_rhs + a43*Yh3_rhs );
                Yh4_rhs = nonlin_h_fun(Yh4, xi, filtr);         
                % 5th stage
                Yh5 = U51*Yh + dt*( a51*Yh1_rhs + a52*Yh2_rhs + ...
                                    a53*Yh3_rhs + a54*Yh4_rhs );
                Yh5_rhs = nonlin_h_fun(Yh5, xi, filtr);                
                % 6th stage
                Yh6 = U61*Yh + dt*( a61*Yh1_rhs + a62*Yh2_rhs + ...
                                    a63*Yh3_rhs + a64*Yh4_rhs + a65*Yh5_rhs );
                Yh6_rhs = nonlin_h_fun(Yh6, xi, filtr);                
                % Solution update
                Yh_new = v11*Yh + dt*( ...
                          b11*Yh1_rhs + ...
                          b13*Yh3_rhs + ...
                          b14*Yh4_rhs + ...
                          b15*Yh5_rhs + ...
                          b16*Yh6_rhs );
            case 'Krogstad' % Krogstad method of order 4 and stiff order 3
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi02*Yh + dt*(1/2)*Phi12*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                Yh3 = Phi02*Yh + dt*( ( (1/2)*Phi12-Phi22 )*Yh1_rhs + Phi22*Yh2_rhs );
                Yh3_rhs = nonlin_h_fun(Yh3,xi,filtr); % function eval in frequency space
                % 4th stage
                Yh4 = Phi0*Yh + dt*( (Phi1-2*Phi2)*Yh1_rhs + 2*Phi2*Yh3_rhs );
                Yh4_rhs = nonlin_h_fun(Yh4,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (Phi1-3*Phi2+4*Phi3)*Yh1_rhs + ...
                                (2*Phi2-4*Phi3)*Yh2_rhs +...
                                (2*Phi2-4*Phi3)*Yh3_rhs +...
                                (-Phi2+4*Phi3)*Yh4_rhs );
            case 'HochbruckOstermann' % 5-stage, order 4, stiff order 4
                % Precompute specific HO coefficients
                % a52 = 1/2*phi2,2 - phi3 + 1/4*phi2 - 1/2*phi3,2
                % a54 = 1/4*phi2,2 - a52
                A52 = 0.5*Phi22 - Phi3 + 0.25*Phi2 - 0.5*Phi32;
                A54 = 0.25*Phi22 - A52;
    
                % 1st stage
                % Yh1 = Yh; % in frequency space
                Yh1_rhs = nonlin_h_fun(Yh,xi,filtr); % function eval in frequency space
                % 2nd stage
                Yh2 = Phi02*Yh + dt*(1/2)*Phi12*Yh1_rhs;
                Yh2_rhs = nonlin_h_fun(Yh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                Yh3 = Phi02*Yh + dt*( ( (1/2)*Phi12-Phi22 )*Yh1_rhs + Phi22*Yh2_rhs );
                Yh3_rhs = nonlin_h_fun(Yh3,xi,filtr); % function eval in frequency space
                % 4th stage
                Yh4 = Phi0*Yh + dt*( (Phi1-2*Phi2)*Yh1_rhs + Phi2*Yh2_rhs + Phi2*Yh3_rhs );
                Yh4_rhs = nonlin_h_fun(Yh4,xi,filtr); % function eval in frequency space
                % 5th stage
                Yh5 = Phi02*Yh + dt*( (0.5*Phi12-2*A52-A54)*Yh1_rhs + A52*Yh2_rhs + A52*Yh3_rhs + A54*Yh4_rhs );
                Yh5_rhs = nonlin_h_fun(Yh5,xi,filtr); % function eval in frequency space
                % Solution update
                Yh_new = Phi0*Yh + dt*( (Phi1 - 3*Phi2 + 4*Phi3)*Yh1_rhs + ...
                                       (-Phi2 + 4*Phi3)*Yh4_rhs + ...
                                       (4*Phi2 - 8*Phi3)*Yh5_rhs );
        end
        Y_new = real(ifft(reshape(Yh_new, 3, N),[],2)); % This applies inverse fft row-wise
        Y = Y_new; tn = tn+dt;
    end
    u = Y(1,:); v = Y(2,:); w = Y(3,:);
    ExpIntTime = toc;
end

%--------------------------------functions--------------------------------%
% Nonlinear part of the RHS
function Yh_rhs = nonlin_h_fun(Yh,xi,filtr)
    Yh_rhs = zeros(length(Yh),1);  % Initialize Y with zeros
    uhat = Yh(1:3:end,1); 
    uh_rhs = -(1j*xi'/2).*fft(real(ifft(uhat.*filtr')).^2);
    Yh_rhs(1:3:end) = uh_rhs;
end
% Frobenius covariants of a given matrix A of size 3 by 3
function [la,mu,nu,Ps,Pp,Pn] = Frob_covs(A)
    EV = eig(A); I = eye(3);
    la = EV(1); mu = EV(2); nu = EV(3);
    Ps = (A-mu*I)*(A-nu*I)/((la-mu)*(la-nu));
    Pp = (A-la*I)*(A-nu*I)/((mu-la)*(mu-nu));
    Pn = (A-la*I)*(A-mu*I)/((nu-la)*(nu-mu));
end
% Computation of phi functions
function [phi_0, phi_1, phi_2, phi_3] = phi_funs_on_mat(dt,L)
    [la,mu,nu,Ps,Pp,Pn] = Frob_covs(L);
    % phi_0
    phi_0 = phi0_z(dt*la)*Ps + phi0_z(dt*mu)*Pp + phi0_z(dt*nu)*Pn;
    % phi_1
    phi_1 = phi1_z(dt*la)*Ps + phi1_z(dt*mu)*Pp + phi1_z(dt*nu)*Pn;
    % phi_2
    phi_2 = phi2_z(dt*la)*Ps + phi2_z(dt*mu)*Pp + phi2_z(dt*nu)*Pn;
    % phi_3
    phi_3 = phi3_z(dt*la)*Ps + phi3_z(dt*mu)*Pp + phi3_z(dt*nu)*Pn;
end
%-----------------------%
% % Phi functions on z
% function u = phi0_z(z)
%     u = exp(z);
% end
% function u = phi1_z(z)
%     u = (exp(z) - 1)./z;
% end
% function u = phi2_z(z)
%     u = (exp(z) - 1 - z)./z.^2;
% end
% function u = phi3_z(z)
%     u = (exp(z) - 1 - z - 0.5*z.^2)./z.^3;
% end
%-----------------------%


% Phi functions on z
function u = phi0_z(z)
    u = exp(z);
end

function [u] = phi1_z(z)
    % Define the threshold for switching (based on machine precision)
    EPS = 1e-8;

    if abs(z) < EPS
        % Case 1: Small |z|, use short Taylor series (k=0 to k=4 terms)
        % Taylor: 1 + z/2! + z^2/3! + z^3/4! + z^4/5!
        u = 1 + z/2 + z^2/6 + z^3/24 + z^4/120;
    else
        % Case 2: Large |z|, use stable formula leveraging expm1(z)
        u = expm1(z) / z;
    end
end

function [u] = phi2_z(z)
    % Threshold for switching to Taylor series (e.g., 1e-8)
    EPS = 1e-8;

    if abs(z) < EPS
        % Case 1: Small |z|, use short Taylor series (k=0 to k=4 terms)
        % Taylor: 1/2! + z/3! + z^2/4! + z^3/5! + z^4/6!
        u = 1/2 + z/6 + z^2/24 + z^3/120 + z^4/720;
    else
        % Case 2: Large |z|, use stable formula leveraging expm1(z)
        % expm1(z) is highly accurate for (e^z - 1)
        u = (expm1(z) - z) / z^2;
    end
end

function [u] = phi3_z(z)
    % Threshold for switching to Taylor series (e.g., 1e-8)
    EPS = 1e-8;
    
    if abs(z) < EPS
        % Case 1: Small |z|, use short Taylor series (k=0 to k=4 terms)
        % Taylor: 1/3! + z/4! + z^2/5! + z^3/6! + z^4/7!
        u = 1/6 + z/24 + z^2/120 + z^3/720 + z^4/5040;
    else
        % Case 2: Large |z|, use stable formula leveraging expm1(z)
        % expm1(z) is highly accurate for (e^z - 1)
        u = (expm1(z) - z - 0.5*z^2) / z^3;
    end
end

%--------------------------------End--------------------------------------%
