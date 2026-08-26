function [dt, u] = CompKdVSolByExpInt(Len, exp_int, tf, N, u0, dt)
%====================================================================%
%  CompKdVSolByExpInt: Solves the KdV equation using exponenetial    %
% time integrators.                                                  %
%--------------------------------------------------------------------%
%  Inputs:                                                           %
%    Len      - length of the computational domain                   %
%    exp_int  - name of the exponential integrator                   %
%    tf       - Final simulation time                                %
%    N        - Number of spatial grid points                        %
%    u0       - Initial condition                                    %
%    dt       - time step                                            %
%                                                                    %
%  Outputs:                                                          %
%    u   - Solution at final time                                    %
% Note: This codes uses two external matlab files: 1. wantcache.m &  %
% phipade.m                                                          %
%====================================================================%    
    % Wave numbers
    xi = -N/2:N/2-1;        % Indices of the wave numbers
    xi = fftshift(xi);      % Reordering of the indices
    xi = xi*2*pi/Len;         % Scaling by 2π/L

    % Filtering
    filtr = ones(size(xi));  % Create an array of ones with the same size as xi
    xi_max = max(abs(xi));   % Compute the maximum absolute value of xi
    filtr(abs(xi) > xi_max * (2/3)) = 0;  % Set elements to 0 where condition is met

    % The initial profile, obtained from analytical solution
    u = u0;
    % Time step 
    nt = round(tf/dt);
    dt = tf/nt; tn = 0;

    % Computing the matrix exponential using scaling and squaring invariants 
    % and Padé approximation. 
    Z = dt*(-(1j*xi).^3);
    % Compute phi0
    phi0 = exp(Z); phi0(1) = 1;
    % Compute phi02
    phi02 = exp(Z/2); phi02(1) = 1;
    % Computing exp using diagonal Pade approximation
    [phi_1, phi_2, phi_3] = phipade(Z, 3);

    % Compute phi1
    phi1 = spdiags(phi_1,0).';
    % Compute phi2
    phi2 = spdiags(phi_2,0).';
    % Compute phi3
    phi3 = spdiags(phi_3,0).';
    % Compute phi12
    [phi_12,phi_22,phi_32] = phipade(Z/2, 3);
    phi12 = spdiags(phi_12,0).';
    phi22 = spdiags(phi_22,0).';
    phi32 = spdiags(phi_32,0).';

    % Compute phi2_hat = phi2_3/5
    [~,phi_2_hat] = phipade(3*Z/5, 2);
    phi2hat = spdiags(phi_2_hat,0).';

    % Precompute ETD5RKF update matrices outside the loop to save time
    if strcmp(exp_int, 'ETD5RKF')
        b11 = (47/150)*phi1 - (188/75)*phi2 + (94/15)*phi3;
        % b12 = zeros
        b13 = (-43/25)*phi1 + (132/5)*phi2 - 66*phi3;
        b14 = (4124/75)*phi1 - (6152/15)*phi2 + (2704/3)*phi3;
        b15 = (189/10)*phi1 - (662/5)*phi2 + 284*phi3;
        b16 = (-1787/25)*phi1 + (12966/25)*phi2 - (5628/5)*phi3;

        %
        c = [0, 2/9, 1/3, 3/4, 1, 5/6];
        u11 = exp(c(1)*Z); u21 = exp(c(2)*Z); u31 = exp(c(3)*Z);
        u41 = exp(c(4)*Z); u51 = exp(c(5)*Z); u61 = exp(c(6)*Z);

        % 
        v11 = phi0;
        
        % Stage coefficients for A(z)
        a21 = -2/3*phi2 + (10/9)*phi2hat;
        a31 = (569/11544)*phi2 + (1355/11544)*phi2hat;
        a32 = (-831/3848)*phi2 + (2755/3848)*phi2hat;
       
        a41 = -(77157/61568)*phi2 + (143535/61568)*phi2hat;
        a42 = (587979/61568)*phi2 - (821745/61568)*phi2hat;
        a43 = (-405/64)*phi2 + (675/64)*phi2hat;
        
        a51 = (655263/7696)*phi2 - (2031205/23088)*phi2hat;
        a52 = (-1148769/7696)*phi2 + (1252665/7696)*phi2hat;
        a53 = (1593/40)*phi2 - (405/8)*phi2hat;
        a54 = (144/5)*phi2 - (80/3)*phi2hat;

        a61 = -(2212835/277056)*phi2 + (6888625/831168)*phi2hat;
        a62 = (477285/30784)*phi2 - (496525/30784)*phi2hat;
        a63 = (-39/16)*phi2 + (65/16)*phi2hat;
        a64 = (-4/9)*phi2 + (20/27)*phi2hat;
        a65 = (-185/96)*phi2 + (575/288)*phi2hat;
    end

    % Time stepping
    for n = 1:nt
        uh = fft(u); 
        % Different exponential tme integrators
        switch exp_int
            %% Lawson methods
            case 'Lawson-Euler' % Lawson-Euler method of order 1 and stiff order 1
                uh_rhs = nonlin_h_fun(uh,xi,filtr);
                uh_new = phi0.*uh + dt*phi0.*uh_rhs;
            case 'Lawson2b' % Lawson method of order 2 and stiff order 1
                % 1st stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi0.*uh + dt*phi0.*uh1_rhs;
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (1/2)*phi0.*uh1_rhs + (1/2)*uh2_rhs );            
            case 'Lawson4' % Lawson method of order 4 and stiff order 1
                % first stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi02.*uh + (dt/2)*phi02.*uh1_rhs;
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                uh3 = phi02.*uh + (dt/2)*uh2_rhs;
                uh3_rhs = nonlin_h_fun(uh3,xi,filtr); % function eval in frequency space
                % 4th stage
                uh4 = phi0.*uh + dt*phi02.*uh3_rhs; 
                uh4_rhs = nonlin_h_fun(uh4,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (1/6)*phi0.*uh1_rhs +...
                    (1/3)*phi02.*uh2_rhs + (1/3)*phi02.*uh3_rhs + ...
                    (1/6)*uh4_rhs );
            %% ETD (Exponential Time Differencing) methods
            case 'Norsett-Euler' % Norsett-Euler method of order 1 and stiff order 1
                uh_rhs = nonlin_h_fun(uh,xi,filtr);
                uh_new = phi0.*uh + dt*phi1.*uh_rhs;  
            case 'ETD2RK' % ETD2RK method of order 2 and stiff order 2
                % first stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi0.*uh + dt*phi1.*uh1_rhs; % in frequency space
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (phi1-phi2).*uh1_rhs + phi2.*uh2_rhs );
            case 'ETD3RK' % ETD2RK method of order 3 and stiff order 2
                % first stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi02.*uh + dt*(1/2)*phi12.*uh1_rhs; % in frequency space
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                uh3 = phi0.*uh + dt*( -phi1.*uh1_rhs + 2*phi1.*uh2_rhs ); % in frequency space
                uh3_rhs = nonlin_h_fun(uh3,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (phi1-3*phi2+4*phi3).*uh1_rhs + (4*phi2-8*phi3).*uh2_rhs + (-phi2+4*phi3).*uh3_rhs );
            case 'ETD4RK' % ETD4RK method of order 4 and stiff order 2
                % first stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi02.*uh + dt*(1/2)*phi12.*uh1_rhs; % in frequency space
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                uh3 = phi02.*uh + dt*(1/2)*phi12.*uh2_rhs; % in frequency space
                uh3_rhs = nonlin_h_fun(uh3,xi,filtr); % function eval in frequency space
                % 4th stage
                uh4 = phi0.*uh + dt*( (1/2)*phi12.*(phi02-1).*uh1_rhs + phi12.*uh3_rhs ); % in frequency space
                uh4_rhs = nonlin_h_fun(uh4,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (phi1-3*phi2+4*phi3).*uh1_rhs +...
                    (2*phi2-4*phi3).*uh2_rhs +...
                    (2*phi2-4*phi3).*uh3_rhs +...
                    (-phi2+4*phi3).*uh4_rhs );
            case 'ETD5RKF'
                % 1st stage
                uh1 = u11.*uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh1, xi, filtr);            
                % 2nd stage
                uh2 = u21.*uh + dt*( a21.*uh1_rhs );
                uh2_rhs = nonlin_h_fun(uh2, xi, filtr);                
                % 3rd stage
                uh3 = u31.*uh + dt*( a31.*uh1_rhs + a32.*uh2_rhs );
                uh3_rhs = nonlin_h_fun(uh3, xi, filtr);                
                % 4th stage
                uh4 = u41.*uh + dt*( a41.*uh1_rhs + a42.*uh2_rhs + a43.*uh3_rhs );
                uh4_rhs = nonlin_h_fun(uh4, xi, filtr);         
                % 5th stage
                uh5 = u51.*uh + dt*( a51.*uh1_rhs + a52.*uh2_rhs + ...
                                    a53.*uh3_rhs + a54.*uh4_rhs );
                uh5_rhs = nonlin_h_fun(uh5, xi, filtr);                
                % 6th stage
                uh6 = u61.*uh + dt*( a61.*uh1_rhs + a62.*uh2_rhs + ...
                                    a63.*uh3_rhs + a64.*uh4_rhs + a65.*uh5_rhs );
                uh6_rhs = nonlin_h_fun(uh6, xi, filtr);                
                % Solution update
                uh_new = v11.*uh + dt*( ...
                          b11.*uh1_rhs + ...
                          b13.*uh3_rhs + ...
                          b14.*uh4_rhs + ...
                          b15.*uh5_rhs + ...
                          b16.*uh6_rhs );
            case 'Krogstad' % 'Krogstad' method of order 4 and stiff order 3
                % first stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi02.*uh + dt*(1/2)*phi12.*uh1_rhs; % in frequency space
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                uh3 = phi02.*uh + dt*( ((1/2)*phi12-phi22).*uh1_rhs + phi22.*uh2_rhs ); % in frequency space
                uh3_rhs = nonlin_h_fun(uh3,xi,filtr); % function eval in frequency space
                % 4th stage
                uh4 = phi0.*uh + dt*( (phi1-2*phi2).*uh1_rhs + 2*phi2.*uh3_rhs ); % in frequency space
                uh4_rhs = nonlin_h_fun(uh4,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (phi1-3*phi2+4*phi3).*uh1_rhs +...
                    (2*phi2-4*phi3).*uh2_rhs +...
                    (2*phi2-4*phi3).*uh3_rhs +...
                    (-phi2+4*phi3).*uh4_rhs );
            case 'HochbruckOstermann' % 5-stage, order 4, stiff order 4
                % Precompute specific HO coefficients
                % a52 = 1/2*phi2,2 - phi3 + 1/4*phi2 - 1/2*phi3,2
                % a54 = 1/4*phi2,2 - a52
                A52 = 0.5*phi22 - phi3 + 0.25*phi2 - 0.5*phi32;
                A54 = 0.25*phi22 - A52;

                % 1st stage
                % uh1 = uh; % in frequency space
                uh1_rhs = nonlin_h_fun(uh,xi,filtr); % function eval in frequency space
                % 2nd stage
                uh2 = phi02.*uh + dt*(1/2)*phi12.*uh1_rhs;
                uh2_rhs = nonlin_h_fun(uh2,xi,filtr); % function eval in frequency space
                % 3rd stage
                uh3 = phi02.*uh + dt*( ( (1/2)*phi12-phi22 ).*uh1_rhs + phi22.*uh2_rhs );
                uh3_rhs = nonlin_h_fun(uh3,xi,filtr); % function eval in frequency space
                % 4th stage
                uh4 = phi0.*uh + dt*( (phi1-2*phi2).*uh1_rhs + phi2.*uh2_rhs + phi2.*uh3_rhs );
                uh4_rhs = nonlin_h_fun(uh4,xi,filtr); % function eval in frequency space
                % 5th stage
                uh5 = phi02.*uh + dt*( (0.5*phi12-2*A52-A54).*uh1_rhs + A52.*uh2_rhs + A52.*uh3_rhs + A54.*uh4_rhs );
                uh5_rhs = nonlin_h_fun(uh5,xi,filtr); % function eval in frequency space
                % Solution update
                uh_new = phi0.*uh + dt*( (phi1 - 3*phi2 + 4*phi3).*uh1_rhs + ...
                                       (-phi2 + 4*phi3).*uh4_rhs + ...
                                       (4*phi2 - 8*phi3).*uh5_rhs );
        end
        u_new = real(ifft(uh_new.*filtr)); % This applies inverse fft rowwise.
        u = u_new; tn = tn + dt; 
    end
end
%---------------------------functions-------------------------------------%
% Nonlinear part of the RHS
function uh_rhs = nonlin_h_fun(uhat,xi,filtr)
    uh_rhs = -(1j*xi/2).*fft(real(ifft(uhat.*filtr)).^2);
end