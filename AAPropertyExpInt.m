%% ===================== PARAMETERS =====================
xL = -40; xR = 40; c = 1.2; tf = 5; 
Tau = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]; 
%Tau = [1e-2, 1e-3]; 

% -------- Choose ONE method here --------
 Mthd_name = 'Lawson2b';
% Mthd_name = 'Lawson4';
% Mthd_name = 'ETD2RK';
% Mthd_name = 'ETD3RK';
% Mthd_name = 'ETD4RK';
% Mthd_name = 'Krogstad';
% Mthd_name = 'HochbruckOstermann';

% -------- Method-specific parameters --------
switch strtrim(Mthd_name)
    case 'Lawson2b'
        N = 2^9; dt_vec = 10.^linspace(-1.6, -3.8, 12); P = 2;
    case 'Lawson4'
        N = 2^9; dt_vec = 10.^linspace(-1.2, -2.5, 10); P = 4;
    case 'ETD2RK'
        N = 2^9; dt_vec = 10.^linspace(-1, -3.8, 12); P = 2;
    case 'ETD3RK'
        N = 2^9; dt_vec = 10.^linspace(-0.5, -3, 12); P = 3;
    case 'ETD4RK'
        N = 2^9; dt_vec = 10.^linspace(-0.5, -2.5, 12); P = 4;
    case 'Krogstad'
        N = 2^9; dt_vec = 10.^linspace(-0.5, -2, 12); P = 4;
    case 'HochbruckOstermann'
        N = 2^9; dt_vec = 10.^linspace(-0.5, -2, 12); P = 4;
end

Len = xR-xL; dx = Len/N; x = xL:dx:xR-dx;

%% ===================== COMPUTE ERRORS =====================
Err = cell(length(Tau),1);
dts = cell(length(Tau),1);

for i = 1:length(Tau)

    filename = sprintf('KdVHRefSolbyPetviashvili/RefSolPetviashvili_c_%.1f_tau_%.1e.mat', c, Tau(i));
    load(filename);

    for k = 1:length(dt_vec)

        fprintf('Method: %s, tau = %.1e, dt = %.2e\n', ...
            Mthd_name, Tau(i), dt_vec(k));

        [u0,v0,w0] = petviashvili_ref_sol(xL,xR,uP,vP,wP,xP,c,x,0);

        [~,dt,u,v,w] = CompKdVHSolByExpInt( ...
            Len, Tau(i), Mthd_name, tf, N, u0, v0, w0, dt_vec(k));

        [u_ref,v_ref,w_ref] = petviashvili_ref_sol(xL,xR,uP,vP,wP,xP,c,x,tf);

        Err{i}(1,k) = norm(u-u_ref)/sqrt(N);
        Err{i}(2,k) = norm(v-v_ref)/sqrt(N);
        Err{i}(3,k) = norm(w-w_ref)/sqrt(N);
        dts{i}(k) = dt;
    end
end

%% ===================== PLOT =====================
if ~exist('Figures', 'dir')
   mkdir('Figures')
end

fig = figure('Name', sprintf('%s Convergence', Mthd_name), ...
       'Color', 'w', ...
       'Units', 'normalized', ...
       'Position', [0.05 0.3 0.6 0.35]);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for comp = 1:3
    nexttile; hold on;

    colors = lines(length(Tau));
    all_err = [];

    % ---- Plot curves ----
    for i = 1:length(Tau)
        err_vals = Err{i}(comp,:);
        all_err = [all_err, err_vals];

        plot(dts{i}, err_vals, '-o', ...
            'Color', colors(i,:), ...
            'MarkerFaceColor', colors(i,:), ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('\\tau = 10^{%d}', round(log10(Tau(i)))));
    end
    % ---- Tight limits ----
    ylim([0.8*min(all_err), 1.2*max(all_err)]);

    % ---- Reference slope ----
    ref_dt = dts{1};
    n = length(ref_dt);

    k = min(5, n);
    mid = ceil(n/2);
    idx_start = max(1, mid - floor(k/2));
    idx_end   = min(n, idx_start + k - 1);
    idx_start = max(1, idx_end - k + 1);

    ref_dt_mid = ref_dt(idx_start:idx_end);

    ref_line = 0.2*(ref_dt_mid.^P) * ...
        (Err{1}(comp,idx_end) / ref_dt(idx_end)^P);

    plot(ref_dt_mid, ref_line, ':k', ...
        'LineWidth', 3, ...
        'DisplayName', sprintf('O(\\Deltat^{%d})', P));

    set(gca, 'XScale','log','YScale','log', ...
        'FontSize',15,'LineWidth',1.2);

    grid minor
    xlabel('\Delta t','FontSize',16);

    if comp == 1
        ylabel('Error','FontSize',16);
        legend('Location','best','FontSize',13);
        title('Convergence for $u$','Interpreter','latex');
    elseif comp == 2
        title('Convergence for $v$','Interpreter','latex');
        set(gca,'YTickLabel',[]);
    else
        title('Convergence for $w$','Interpreter','latex');
        set(gca,'YTickLabel',[]);
    end

    box on;
end
% Save figures
filename = sprintf('Figures/%s_AA_c_%.1f_tf_%d.pdf', ...
    Mthd_name, c, tf);
exportgraphics(fig, filename, 'ContentType','vector');

%% ===================== Reference solution =====================
function [ut,vt,wt] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, t)
    xt = mod(x - c*t - xL, xR - xL) + xL;
    ut = spline(xP, uP, xt); 
    vt = spline(xP, vP, xt);
    wt = spline(xP, wP, xt);
end