clc
clear all
close all

% Adding path to utils functions
addpath('../utils')

% Creating path for results
res_path = "1";
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

% Initial condition
h0 = @(x) 1 + 0.5*sin(pi*x);

% Initial data set
IC = @(x) [h0(x); u * h0(x)];
source = @(xc,time,u,g,h0) source_function(xc,time,u,g);
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'periodic';
limiters = ["NONE" "MUSCL" "MINMOD" "TVB"];

% Define different values of \delta x
dx_values = [0.001 0.005 0.0075 0.01];
err_lf = [];
err_roe = [];
dx_iter = 0;

for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    dx_iter = dx_iter + 1;
    
    err_dx_lf = [];
    err_dx_roe = [];
    xc = a+0.5*dx:dx:b-0.5*dx;
    q_exact = find_exact(h0, u, xc, Tfinal);
    q_lf = {};
    q_roe = {};
    
    for lim = limiters
        [xc, q_lf{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, lf_flux, lim, source);
        [xc, q_roe{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, roe_flux, lim, source);
    
        err_dx_lf(end+1) = norm(q_exact - q_lf{end},2);
        err_dx_roe(end+1) = norm(q_exact - q_roe{end},2);
    end
    
    % LF flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 1, xc, q_exact, limiters, q_lf, "Lax-Freidrich", res_path, 0.4, 1.6, 0, 0.5);
    % Roe flux dx plot
    plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 2, xc, q_exact, limiters, q_roe, "Roe", res_path, 0.4, 1.6, 0, 0.5);
    
    err_lf = [err_lf; err_dx_lf];
    err_roe = [err_roe; err_dx_roe];
end

% LF error
dx_iter = dx_iter + 1;
fig = figure(2*(dx_iter - 1) + 1);
% axes('XScale', 'log', 'YScale', 'log')
for i = 1:length(limiters)
    hold all
    txt = ['Limiter ',num2str(limiters(i))];
    loglog(dx_values,err_lf(:, i), 'LineWidth',2, 'DisplayName',txt);
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^3'];
loglog(dx_values,dx_values.^3, 'k-', 'LineWidth',1, 'DisplayName',txt);
legend('Location','NorthWest')
title("Lax-Freidrich flux")
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "LF_error.png");

% Roe error
fig = figure(2*(dx_iter - 1) + 2);
% axes('XScale', 'log', 'YScale', 'log')
for i = 1:length(limiters)
    hold all
    txt = ['Limiter ',num2str(limiters(i))];
    loglog(dx_values,err_roe(:, i), 'LineWidth',2, 'DisplayName',txt); 
end
txt = ['h'];
loglog(dx_values,dx_values, 'k-.', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^3'];
loglog(dx_values,dx_values.^3, 'k-', 'LineWidth',1, 'DisplayName',txt);
title("Roe flux")
legend('Location','NorthWest')
grid on;
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "roe_error.png");