maxNumCompThreads(1);
xL = -40; xR = 40; c = 1.2; tf = 1; tau = 1e-5;
% Choose an order
p = 4;
switch p 
    case 2
        Exp_Mthds = {'ETD2RK'};
        ImEx_Mthds = {'AGSA(3,4,2)'}; Type = {'I'};
        N = 2^10;
        dts = 10.^linspace(-1.5,-3.5,10);
        Len = xR-xL; dx = Len/N; x = xL:dx:xR-dx;
    case 3
        Exp_Mthds = {'ETD3RK'};
        ImEx_Mthds = {'ARK3(2)4L[2]SA','ARS(4,4,3)'}; Type = {'II','II'};
        N = 2^10;
        dts = 10.^linspace(-0.5,-3,10);
        Len = xR-xL; dx = Len/N; x = xL:dx:xR-dx;
    case 4
        % 4th order
        Exp_Mthds = {'ETD4RK','HochbruckOstermann'};
        ImEx_Mthds = {'ARK4(3)6L[2]SA','ARK4(3)7L[2]SA'}; Type = {'II','II'};
        N = 2^10;
        dts = 10.^linspace(-1, -2, 10); 
        Len = xR-xL; dx = Len/N; x = xL:dx:xR-dx;
end

%----------Compute------------%
Exp_Err_u = cell(1,length(Exp_Mthds)); ImEx_Err_u = cell(1,length(ImEx_Mthds));
Exp_dts = cell(1,length(Exp_Mthds)); ImEx_dts = cell(1,length(ImEx_Mthds));
Exp_RT = cell(1,length(Exp_Mthds)); ImEx_RT = cell(1,length(ImEx_Mthds));

% load the reference solution, computed on N = 2^10 grid points
% using Petviashvili method with a tol = 1e-12
filename = sprintf('KdVHRefSolbyPetviashvili/RefSolPetviashvili_c_%.1f_tau_%.1e.mat', c, tau);
load(filename);

% Number of runs for filtering CPU timing noise
num_runs = 5; 

% Exponential methods
for j = 1:length(Exp_Mthds)
    for k = 1:length(dts)
        fprintf('Method: %s, τ = %1.1e, and  dt = %.e.\n', Exp_Mthds{j}, tau, dts(k))
        [u0,v0,w0] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, 0);

        % Time-averaging loop to eliminate zig-zags
        run_times = zeros(1, num_runs);
        for r = 1:num_runs
            [RunTime,~,u,v,w] = CompKdVHSolByExpInt(Len, tau, Exp_Mthds{j}, tf, N, u0, v0, w0, dts(k));
            run_times(r) = RunTime;
        end

        [u_ref,v_ref,w_ref] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, tf);
        err_u = norm(u-u_ref)/sqrt(N); 

        Exp_Err_u{1,j}(k) = err_u;
        Exp_dts{1,j}(k) = dts(k);
        Exp_RT{1,j}(k) = mean(run_times, 'omitnan');
    end
end 

% ImEx methods
for j = 1:length(ImEx_Mthds)
    for k = 1:length(dts)
        fprintf('Method: type %s %s, τ = %.1e, and  dt = %.e.\n', Type{j}, ImEx_Mthds{j}, tau, dts(k))
        [u0,v0,w0] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, 0);

        % Time-averaging loop to eliminate zig-zags
        run_times = zeros(1, num_runs);
        for r = 1:num_runs
            [RunTime,~,u,v,w] = CompKdVHSolByImExInt(Len, tau, ImEx_Mthds{j}, Type{j}, tf, N, u0, v0, w0, dts(k));
            run_times(r) = RunTime;
        end

        [u_ref,v_ref,w_ref] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, tf);
        err_u = norm(u-u_ref)/sqrt(N); 

        ImEx_Err_u{1,j}(k) = err_u;
        ImEx_dts{1,j}(k) = dts(k);
        ImEx_RT{1,j}(k) = mean(run_times, 'omitnan');
    end
end
%-------------------------------%
if ~exist('Figures', 'dir')
   mkdir('Figures')
end

%plot
C = {'b','r','m','k','g'}; Cref = {[0.5,0.5,0.5]}; sol_comp = {'u','v','w'};
linS = {':','--','-','-.'}; Mar = {'o','s','*','+','d'}; ms = 10; fs = 20;

%---------for minimum white space
% Ensure figure remains valid
fig = figure(1);
% Set figure size
set(fig, 'Units', 'pixels', 'Position', [50, 50, 500, 500]);
% Adjust axis to remove extra margins
ax = gca;
set(ax, 'LooseInset', [0, 0, 0, 0]);
%---------for minimum white space
legendEntries = {};

for i = 1:length(Exp_Mthds)
    loglog(Exp_RT{1,i}, Exp_Err_u{1,i}, 'linestyle',linS{i},'color',C{i},'marker',Mar{i},'MarkerSize',ms, 'LineWidth', 3)
    legendEntries{end+1} = sprintf('%s', Exp_Mthds{i});
    hold on
end

lexp = length(Exp_Mthds);

for i = 1:length(ImEx_Mthds)
    loglog(ImEx_RT{1,i}, ImEx_Err_u{1,i}, 'linestyle',linS{lexp+i},'color',C{lexp+i},'marker',Mar{lexp+i},'MarkerSize',ms, 'LineWidth', 3)
    legendEntries{end+1} = sprintf('%s', ImEx_Mthds{i});
    hold on
    % Set legend
    legend(legendEntries, 'Interpreter', 'latex','NumColumns',1,'Location','best','Box','off','FontSize',fs)
    % x and ylim
    % xlim()
    % ylim()
    xlabel('Runtime(s)');
    ylabel('Error', 'Rotation', 90);
    title(['$ c = ', sprintf('%.1f', c),', \tau = ',sprintf('\\tau = 10^{%d}', round(log10(tau))),...
        ', t_f = ', sprintf('%d', tf), '$'], 'Interpreter', 'latex');
    set(gca, 'FontSize', 20)
end

% save figures
filename = sprintf('Figures/Pet_Runtime_Exp_vs_ImEx_tau_%.1e_order_%d.pdf',tau,p);  % Create the filename
% Save as PDF: %---------for minimum white space
exportgraphics(gca, filename, 'ContentType', 'vector', 'Resolution', 300);

%---------------------------functions-------------------------------------%
% Reference solution of KdVH by Petviashvili
function [ut,vt,wt] = petviashvili_ref_sol(xL, xR, uP, vP, wP, xP, c, x, t)
    xt = mod(x-c*t-xL,xR-xL)+xL;
    ut = spline(xP,uP,xt); 
    vt = spline(xP,vP,xt);
    wt = spline(xP,wP,xt);
end
%--------------------------------End--------------------------------------%