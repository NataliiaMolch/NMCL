clc
clear all
close all

% Adding path to utils functions
addpath('../utils')

% Creating path for results
res_path = "2_2";
mkdir(res_path);

% Method constants
g = 1;
u = 0.25;
M = 50;

% Discretization and time parameters
Tfinal = 2;
CFL    = 0.5;

% Domain parameters
a = 0;
b = 2;

% Initial data set
IC = @(x) [1 - 0.2*sin(2*pi*x); 0.5*ones(1, length(x))];
source = @(xc,time,u,g,h0) zeros(2, length(xc));
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'periodic';
limiters = ["NONE" "MUSCL" "MINMOD" "TVB"];

% Define different values of \delta x
dx_values = [0.05, 0.01];
err_lf = cell(length(limiters));
err_roe = cell(length(limiters));
dx_iter = 0;

% Reference solution
ref_sol_file = 'ref_sol_2.mat';
load(ref_sol_file,'ref_sol_2_dx')
load(ref_sol_file,'ref_sol_2_x')
load(ref_sol_file,'ref_sol_2')

for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    dx_iter = dx_iter + 1;
    
    q_lf = {};
    q_roe = {};
    
    for lim = 1:length(limiters)
        [xc, q_lf{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, lf_flux, limiters(lim), source);
        [xc, q_roe{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, roe_flux, limiters(lim), source);
        
        % Correlate mesh and fine mesh
        ref_sol = [];
        i = 1; j = 1;
        while i <= length(xc)
            j_start = j;
            while (ref_sol_2_x(j) + 0.5 * ref_sol_2_dx) < i * dx
                j = j + 1;
            end
            ref_sol(:, end + 1) = mean(ref_sol_2(:, j_start:j), 2);
            i = i + 1;
        end
        
        for i=1:2
            err_lf{lim}(i, end+2-i) = norm(ref_sol(i, :) - q_lf{end}(i, :)) / length(q_lf{end}(i, :));
            err_roe{lim}(i, end+2-i) = norm(ref_sol(i, :) - q_roe{end}(i, :)) / length(q_roe{end}(i, :));
        end
    end
    
    % LF flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 1, ref_sol_2_x, ref_sol_2, limiters, q_lf, "Lax-Freidrich", res_path, 0.6, 1.4, 0.2, 0.8);
    % Roe flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 2, ref_sol_2_x, ref_sol_2, limiters, q_roe, "Roe", res_path, 0.6, 1.4, 0.2, 0.8);
end

% LF error
fig = figure();
set(gcf,'position',[10,10,800,400]);
hold on
for i = 1:length(limiters)
    txt = ['Limiter ',num2str(limiters(i))];
    r = polyfit(log(dx_values), log(err_lf{i}(1, :)), 1);
    loglog(dx_values, dx_values.^r(1).*exp(r(2)), 'LineWidth',2, 'DisplayName',txt);
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
legend('Location','NorthEastoutside')
title("Lax-Freidrich flux depth")
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "LF_depth_error.png");

fig = figure();
set(gcf,'position',[10,10,800,400]);
hold on
for i = 1:length(limiters)
    txt = ['Limiter ',num2str(limiters(i))];
    r = polyfit(log(dx_values), log(err_lf{i}(2, :)), 1);
    loglog(dx_values,dx_values.^r(1).*exp(r(2)), 'LineWidth',2, 'DisplayName',txt);
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
legend('Location','NorthEastoutside')
title("Lax-Freidrich flux discharge")
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "LF_discharge_error.png");

% Roe error
fig = figure();
set(gcf,'position',[10,10,800,400]);
hold on
for i = 1:length(limiters)
    txt = ['Limiter ',num2str(limiters(i))];
    r = polyfit(log(dx_values), log(err_roe{i}(1, :)), 1);
    loglog(dx_values,dx_values.^r(1).*exp(r(2)), 'LineWidth',2, 'DisplayName',txt);
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
legend('Location','NorthEastoutside')
title("Roe flux depth")
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "Roe_depth_error.png");

fig = figure();
set(gcf,'position',[10,10,800,400]);
hold on
for i = 1:length(limiters)
    txt = ['Limiter ',num2str(limiters(i))];
    r = polyfit(log(dx_values), log(err_roe{i}(2, :)), 1);
    loglog(dx_values,dx_values.^r(1).*exp(r(2)), 'LineWidth',2, 'DisplayName',txt);
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
legend('Location','NorthEastoutside')
title("Roe flux discharge")
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "Roe_discharge_error.png");