xL = -40; xR = 40; c = 1.2; tf = 0.01; tau = 1e-4; % also use ETD2/3 ETD5RKF
N = 2^9; dt = 10^(-3); Len = xR-xL; dx = Len/N; x = xL:dx:xR-dx; exp_int = 'ETD2RK';

xi = (-N/2:N/2-1);        % Indices of the wave numbers
xi = fftshift(xi);      % Reordering of the indices
xi = xi*2*pi/Len;         % Scaling by 2π/L

filename = sprintf('KdVHRefSolbyPetviashvili/RefSolPetviashvili_c_%.1f_tau_%.1e.mat', c, tau);
load(filename);
[u0,v0,w0] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, 0);

% Figure names
fig_names ={'wellprepared','unprepared','zero_order_wp','first_order_wp'};

% choose option for WP = {'wp','nwp','wp0','wp1'}
WP = 'nwp';

switch WP
    case 'wp'
        % Case 1: well prepared initial data
        yscale = 1e-12;
        tstring = 'Well prepared';
        fprintf('\n %s \n', tstring);
        fig_name = fig_names{1};
    case 'nwp'
        % Case 2: not well prepared initial data
        v0(1:end) = 0.;
        w0(1:end) = 0.;
        yscale = 1.e-3;
        tstring = 'Not well prepared';
        fprintf('\n %s \n', tstring);
        fig_name = fig_names{2};
    case 'wp0'
        % Case 3: zero-order well prepared initial data
        v0 = ifft(1j*xi.*fft(u0));
        w0 = ifft(1j*xi.*fft(v0));
        yscale = 1.e-7;
        tstring = "Prepared to order zero";
        fprintf('\n %s \n', tstring);
        fig_name = fig_names{3};
    case 'wp1'
        % Case 4: first-order well prepared initial data
        u0x = ifft(1j*xi.*fft(u0));
        w0x = ifft((1j*xi).^3.*fft(u0));
        v0t = -ifft(1j*xi.*fft(u0.*u0x+w0x));
        v1 = ifft(1j*xi.*fft(v0t));
        w1 = ifft(1j*xi.*fft(v1))-v0t;
        v0 = ifft(1j*xi.*fft(u0)) + tau*v1;
        w0 = ifft(-xi.^2.*fft(u0)) + tau*w1;
        yscale = 1.e-11;
        tstring = "Prepared to first order";
        fprintf('\n %s \n', tstring);
        fig_name = fig_names{4};
end

[ExpTime,dt,u,v,w] = CompKdVHSolByExpInt(Len, tau, exp_int, tf, N, u0, v0, w0, dt);

% reference solution at the final time
[u_ref,v_ref,w_ref] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP,c,x,tf);

[dt,uSol] = CompKdVSolByExpInt(Len, exp_int, tf, N, u0, dt);

plot(x,u,'-b',x,u_ref,':g',x,uSol,'--k','linewidth',2)
title(tstring)
ylabel('u'); xlabel('x');
ylim([-yscale, yscale])
grid minor
legend('u','u_{exact}','KdV Sol','Location','best')
set(gca,'FontSize',14)

% save figures
if ~exist('Figures', 'dir')
   mkdir('Figures')
end

filename = sprintf('Figures/%s.pdf', fig_name);
exportgraphics(gca, filename, 'ContentType', 'vector', 'Resolution', 300);

%---------------------------functions-------------------------------------%
function [ut,vt,wt] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, t)
    xt = mod(x-c*t-xL,xR-xL)+xL;
    ut = spline(xP,uP,xt); 
    vt = spline(xP,vP,xt);
    wt = spline(xP,wP,xt);
end
%--------------------------------End--------------------------------------%