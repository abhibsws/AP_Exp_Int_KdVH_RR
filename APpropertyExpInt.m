%=========================================================================%
% This code produces numerical results for the Asymptotic Preserving (AP), 
% plots and the generates LaTeX tables with Order of Convergence.
%=========================================================================%
% --- Configuration Parameters ---
xL = -40; xR = 40; 
c = 1.2; 

dts = [0.003,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.0001];
dts_used = zeros(1,length(dts));
% Exponential methods
Mthd_names = {'Lawson-Euler','Lawson2b','Lawson4','Norsett-Euler',...
'ETD2RK','ETD3RK','ETD4RK','Krogstad','HochbruckOstermann','ETD5RKF'};
tf = 5;

% For experiment
% dts = [0.015]; dts_used = zeros(1,length(dts));
% Mthd_names = {'HochbruckOstermann'}; tf = 5;

%dts = [0.0001]; dts_used = zeros(1,length(dts)); 
% Mthd_names = {'ETD5RKF'}; tf = 5;

% dts = [0.015]; dts_used = zeros(1,length(dts));
% Mthd_names = {'ETD4RK'}; tf = 5;

% dts = [0.015]; dts_used = zeros(1,length(dts));
% Mthd_names = {'Krogstad'}; tf = 5;




N = 2^9; % Define N outside the loop for use in plotting titles
Len = xR-xL; dx = Len/N; x = (xL:dx:xR-dx);

% Relaxation parameters (tau is decreasing by factor of 10)
Tau = [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10];

% --- Initialize Error Storage ---
Err_u = zeros(length(Mthd_names),length(Tau)); 
Err_v = zeros(length(Mthd_names),length(Tau));
Err_w = zeros(length(Mthd_names),length(Tau));

% --- Generating Initial and Reference Solutions ---
% generating IC on the grid using ref solution
u0 = initial_profile(c,x); % a row vector
% Well prepared initial data for v and w upto leading order
D = fourier_first_derivative_matrix(xL, xR, N);
v0 = (D*u0')'; w0 = (D^2*u0')'; % row vectors

% --- Main Computation Loop ---
for i = 1:length(Mthd_names)
    % Compute the numerical solution of the kdV equation (Reference solution)
    [~,uSol] = CompKdVSolByExpInt(Len, Mthd_names{i}, tf, N, u0, dts(i));
    u_KdV = uSol;

    % Generating reference solution for v and w component for KdV numerical sol
    v_KdV = (D*uSol')'; w_KdV = (D^2*uSol')';

    for j = 1:length(Tau)
        fprintf('Method: %s and tau = %.1e.\n', Mthd_names{i}, Tau(j))

        % computing solution by an exp integrator for the KdVH system 
        [~,dt,u,v,w] = CompKdVHSolByExpInt(Len, Tau(j), Mthd_names{i}, tf, N, u0, v0, w0, dts(i));

        % Compute error (L2 norm)
        Err_u(i,j) = norm(u-u_KdV)/sqrt(N); 
        Err_v(i,j) = norm(v-v_KdV)/sqrt(N);
        Err_w(i,j) = norm(w-w_KdV)/sqrt(N);
        dts_used(1,i) = dt;
    end
end


% -------------------------------------------------------------------------
% --- LaTeX Table Generation and EOC Calculation ---
% -------------------------------------------------------------------------

for i = 1:length(Mthd_names)
    method_name = Mthd_names{i};

    % 1. Calculate EOC
    EOC_u = NaN(1, length(Tau));
    EOC_v = NaN(1, length(Tau));
    EOC_w = NaN(1, length(Tau));

    % Loop from the second step (j=2) to calculate EOC relative to the previous step (j-1)
    for j = 2:length(Tau)
        % The step size ratio (Tau_j-1 / Tau_j) is 10 for this setup, so log10(ratio) = 1
        tau_ratio_log = log10(Tau(j-1) / Tau(j)); 

        % EOC = log(Error_prev / Error_curr) / log(Tau_prev / Tau_curr)
        EOC_u(j) = log10(Err_u(i, j-1) / Err_u(i, j)) / tau_ratio_log;
        EOC_v(j) = log10(Err_v(i, j-1) / Err_v(i, j)) / tau_ratio_log;
        EOC_w(j) = log10(Err_w(i, j-1) / Err_w(i, j)) / tau_ratio_log;
    end

    % save
    if ~exist('AP_Data', 'dir')
       mkdir('AP_Data')
    end


    % 2. Generate LaTeX file
    safe_method_name = strrep(method_name, '-', '_'); % Replace dashes for filename
    filename = sprintf('AP_Data/AP_data_%s.tex', safe_method_name);

    fileID = fopen(filename, 'w');

    % LaTeX Table Header
    fprintf(fileID, '\\begin{table}[ht]\n');
    fprintf(fileID, '\\centering\n');
    fprintf(fileID, '\\caption{Asymptotic Preserving property of the %s method for the KdVH system. The $\\ell_2$ norms of the errors are calculated relative to the numerical solution of the KdV equation.}\n', method_name);
    fprintf(fileID, '\\label{tab:ap_results_%s}\n', safe_method_name);
    fprintf(fileID, '\\resizebox{\\textwidth}{!}{\\begin{tabular}{ccccccc}\n'); % Use resizebox to ensure fit
    fprintf(fileID, '  \\hline\n');
    fprintf(fileID, '  \\textbf{$\\tau$} & \\textbf{error u} & \\textbf{EOC u} & \\textbf{error v} & \\textbf{EOC v} & \\textbf{error w} & \\textbf{EOC w} \\\\ \\hline\n');

    % 3. Write Data Rows
    for j = 1:length(Tau)
        tau_val = Tau(j);
        err_u = Err_u(i, j);
        err_v = Err_v(i, j);
        err_w = Err_w(i, j);

        % Format EOCs, using NaN check for the first entry
        eoc_u_str = format_value(EOC_u(j), 'EOC');
        eoc_v_str = format_value(EOC_v(j), 'EOC');
        eoc_w_str = format_value(EOC_w(j), 'EOC');

        % Format the row: Tau | Error u | EOC u | Error v | EOC v | Error w | EOC w
        fprintf(fileID, '  %1.2e & %1.2e & %s & %1.2e & %s & %1.2e & %s \\\\\n', ...
            tau_val, err_u, eoc_u_str, err_v, eoc_v_str, err_w, eoc_w_str);
    end

    % LaTeX Table Footer
    fprintf(fileID, '  \\hline\n');
    fprintf(fileID, '\\end{tabular}}\n');
    fprintf(fileID, '\\end{table}\n');

    fclose(fileID);
    fprintf('Successfully generated LaTeX table: %s\n', filename);
end

% -------------------------------------------------------------------------
% --- PLOTTING ---
% -------------------------------------------------------------------------
C = {'b','r','g','k','m'}; Cref = {[0.5,0.5,0.5]};
linS = {':','--','-'}; Mar = {'o','s','+','*','d'}; ms = 10; fs = 22;
ST = [3,3,3,3,3,3,3,3,3,3]; EN = [6,6,6,6,6,6,6,6,6,6]; Coeff = [1e0,1e0,1e0,1e0,1e0,1e0,1e0,1e0,1e0,1e0];
YL = [1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8]; 
YU = [2e-1,2e-1,2e-1,5e-3,5e-3,5e-3,5e-3,5e-3,5e-3,5e-3];

for i = 1:length(Mthd_names)
    %---------for minimum white space
    fig = figure(i);
    set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
    ax = gca;
    set(ax, 'LooseInset', [0, 0, 0, 0]);
    %---------for minimum white space
    
    loglog(Tau, Err_u(i,:), 'o-b', 'LineWidth', 2)
    hold on
    loglog(Tau, Err_v(i,:), 's-r', 'LineWidth', 2)
    hold on
    loglog(Tau, Err_w(i,:), '*-m', 'LineWidth', 2)
    hold on
    loglog(Tau(ST(i):EN(i)), Coeff(i) * Tau(ST(i):EN(i)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
    legend('u','v','w', '$\mathcal{O}(\tau)$', 'FontSize', fs, 'Location', 'best', 'Interpreter', 'latex')
    grid minor
    xlabel_handle = xlabel('$\tau$', 'Interpreter', 'latex');
    xlabel_handle.Position(2) = xlabel_handle.Position(2) + 0.15*YL(i);
    ylabel('Error', 'Interpreter', 'latex')
    xlim([Tau(end), Tau(1)]); 
    % ylim([YL(i), YU(i)]); % Commented out to prevent errors with some versions
    
    title([sprintf('%s: ', Mthd_names{i}), 'c = ', sprintf('%.1f', c),...
        '$, N = ', sprintf('%d', N), ', dt = ', sprintf('%.3f', dts(i)),...
        ', t_f = ', sprintf('%d', tf), '$'], 'Interpreter', 'latex');
    set(gca, 'FontSize', fs-4) 
    
    % save
    if ~exist('Figures', 'dir')
       mkdir('Figures')
    end
    filename = sprintf('Figures/AP_c_%.1f_N_%d_tf_%d_%s.pdf',c,N,tf,Mthd_names{i});
    %exportgraphics(gca, filename, 'ContentType', 'vector', 'Resolution', 300);
end
%---------------------------AUXILIARY FUNCTIONS----------------------------%

% Helper function to format EOC values
function str = format_value(val, type)
    if isnan(val)
        str = 'NaN';
    elseif strcmp(type, 'EOC')
        % Format EOC to two decimal places
        str = sprintf('%1.2f', val);
    else
        % Default scientific notation for errors
        str = sprintf('%1.2e', val);
    end
end

% Initial condition
function u0 = initial_profile(c,x)
    mu = 0; A = 3*c;
    u0 = A*sech(sqrt(3*A)/6*(x-mu)).^2;  
end
% Fourier differetiation matrix
function [D] = fourier_first_derivative_matrix(xmin, xmax, Nx)
    Len = xmax - xmin; 
    x = linspace(xmin, xmax, Nx+1)'; ix = x(1:end-1); Nf = Nx;
    
    % To compute fhat
    omega = zeros(Nf, Nx);
    for ii = 1 : Nf
        kk = ii - Nf / 2;
        for jj = 1 : Nx
            omega(ii, jj) = exp(-1i*kk*ix(jj) * 2 * pi / Len) / Nx;
        end
    end
    
    % Now take the inverse using fhat
    iomega = zeros(Nx, Nf);
    for ii = 1 : Nx
        for jj = 1 : Nf
            kk = jj - Nf / 2;
            iomega(ii, jj) = 1i * kk  * 2 * pi / Len * exp(1i * kk * ix(ii) * 2 * pi / Len);
        end
    end
    
    % Cheating ;-)
    D = real(iomega * omega);
end
%--------------------------------End--------------------------------------%
