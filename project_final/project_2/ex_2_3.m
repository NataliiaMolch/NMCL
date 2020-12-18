clc
clear all
close all

% Adding path to utils functions
addpath('../utils')

% Creating path for results
res_path = "3";
mkdir(res_path);

% Method constants
g = 1;
u = 0.25;
M = 10;

% Discretization and time parameters
Tfinal = 2;
CFL    = 0.5;

% Domain parameters
a = 0;
b = 2;

% Initial data set
IC = @(x) [ones(1, length(x)); (x < 1) * (-1.5)];
source = @(xc,time,u,g,h0) zeros(2, length(xc));
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'open';
limiters = ["NONE" "MUSCL" "MINMOD" "TVB"];

% Define different values of \delta x
dx_values = [0.01, 0.005, 0.001];
err_lf = [];
err_roe = [];
dx_iter = 0;

% Reference solution
ref_sol_file = 'ref_sol_3.mat';
load(ref_sol_file,'ref_sol_3_dx')
load(ref_sol_file,'ref_sol_3_x')
load(ref_sol_file,'ref_sol_3')

for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    dx_iter = dx_iter + 1;
    
    err_dx_lf = [];
    err_dx_roe = [];
    q_lf = {};
    q_roe = {};
    
    for lim = limiters
        [xc, q_lf{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, lf_flux, lim, source);
        [xc, q_roe{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, roe_flux, lim, source);
        
        % Correlate mesh and fine mesh
        ref_sol = [];
        i = 1; j = 1;
        while i <= length(xc)
            j_start = j;
            while (ref_sol_3_x(j) + 0.5 * ref_sol_3_dx) < i * dx
                j = j + 1;
            end
            ref_sol(:, end + 1) = mean(ref_sol_3(:, j_start:j), 2);
            i = i + 1;
        end
        
        err_dx_lf(end+1) = norm(ref_sol - q_lf{end},2);
        err_dx_roe(end+1) = norm(ref_sol - q_roe{end},2);
    end
    
    % LF flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 1, ref_sol_3_x, ref_sol_3, limiters, q_lf, "Lax-Freidrich", res_path, 0.2, 1.2, -1, 0.2);
    % Roe flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 2, ref_sol_3_x, ref_sol_3, limiters, q_roe, "Roe", res_path, 0.2, 1.2, -1, 0.2);
    
    err_lf = [err_lf; err_dx_lf];
    err_roe = [err_roe; err_dx_roe];
end

% LF error
dx_iter = dx_iter + 1;
fig = figure(2*(dx_iter - 1) + 1);
for i = 1:length(limiters)
    hold all
    txt = ['Limiter ',num2str(limiters(i))];
    plot(dx_values,err_lf(:, i), 'LineWidth',2, 'DisplayName',txt);
end
legend('Location','NorthWest')
title("Lax-Freidrich flux")
grid on;
ylabel('Error')
xlabel('dx')
hold off
saveas(fig, res_path + "/" + "LF_error.png");

% Roe error
fig = figure(2*(dx_iter - 1) + 2);
for i = 1:length(limiters)
    hold all
    txt = ['Limiter ',num2str(limiters(i))];
    plot(dx_values,err_roe(:, i), 'LineWidth',2, 'DisplayName',txt);
end
title("Roe flux")
legend('Location','NorthWest')
grid on;
ylabel('Error')
xlabel('dx')
hold off
saveas(fig, res_path + "/" + "roe_error.png");