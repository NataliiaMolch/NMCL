% M tunning for the step 0.05

clc
clear all
close all

% Adding path to utils functions
addpath('../utils')

% Creating path for results
res_path = "M_1";
mkdir(res_path);

% Method constants
g = 1;
u = 0.25;

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
lim = "TVB";
M_values = [1 5 7 10 15 20 40 50 100 200];
dx =0.005;

% Define different values of \delta x
err_lf = [];
err_roe = [];
dx_iter = 0;

for M = M_values
    disp("Computing M = " + num2str(M));
    dx_iter = dx_iter + 1;
    
    err_dx_lf = [];
    err_dx_roe = [];
    xc = a+0.5*dx:dx:b-0.5*dx;
    q_exact = find_exact(h0, u, xc, Tfinal);
    N = length(q_exact);
    q_lf = {};
    q_roe = {};
    
    
    [xc, q_lf{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, lf_flux, lim, source);
    [xc, q_roe{end+1}] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, roe_flux, lim, source);

    err_dx_lf(end+1) = norm(q_exact - q_lf{end},2)/N;
    err_dx_roe(end+1) = norm(q_exact - q_roe{end},2)/N;
    
    
    % LF flux dx plot
%     plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 1, xc, q_exact, [lim], q_lf, "Lax-Freidrich", res_path, 0.4, 1.6, 0, 0.5);
%     % Roe flux dx plot
%     plot_support(a, b, dx, xc, Tfinal, 2*(dx_iter - 1) + 2, xc, q_exact, [lim], q_roe, "Roe", res_path, 0.4, 1.6, 0, 0.5);
    
    err_lf = [err_lf; err_dx_lf];
    err_roe = [err_roe; err_dx_roe];
end

% LF error
dx_iter = dx_iter + 1;
fig = figure(2*(dx_iter - 1) + 1);
% axes('XScale', 'log', 'YScale', 'log')
hold all
err_lf = err_lf';
plot(M_values,err_lf, 'LineWidth',2);
title("Lax-Freidrich flux")
grid on;
ylabel('Error')
xlabel('M')
hold off
saveas(fig, res_path + "/" + "LF_error.png");

% Roe error
fig = figure(2*(dx_iter - 1) + 2);
% axes('XScale', 'log', 'YScale', 'log')
hold all
err_roe = err_roe';
plot(M_values,err_roe, 'LineWidth',2); 
title("Roe flux")
grid on;
ylabel('Error')
xlabel('M')
hold off
saveas(fig, res_path + "/" + "roe_error.png");