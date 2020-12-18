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
k = 2;

% Discretization and time parameters
Tfinal = 2;
CFL    = 0.5;

% Domain parameters
a = 0;
b = 2;

% Initial condition
h0 = @(x) 1 + 0.5*sin(pi*x);

% Initial data set
ICh = @(x) h0(x); 
ICm = @(x) u * h0(x);
integrated_source = @(xf,time,u,g) source_cell_average(xf,time,u,g);
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'periodic';

% Define different values of \delta x
dx_values = [0.08, 0.006, 0.005, 0.004, 0.002];
err_lf = [];
err_roe = [];
dx_iter = 0;

for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    dx_iter = dx_iter + 1;
    
    xc = a+0.5*dx:dx:b-0.5*dx;
    q_exact = find_exact(h0, u, xc, Tfinal);

    [xc, q_lf] = SSPRK3(a, b, dx, bc, ICh, ICm, u, g, k, Tfinal, CFL, lf_flux, integrated_source);
    [xc, q_roe] = SSPRK3(a, b, dx, bc, ICh, ICm, u, g, k, Tfinal, CFL, roe_flux, integrated_source);
    
    err_lf(end+1) = norm(q_exact - q_lf,2);
    err_roe(end+1) = norm(q_exact - q_roe,2);
    
    fig = figure(dx_iter);
    subplot(2,1,1)
    plot(xc,q_lf(1,:),'-r','LineWidth',2);
    hold all
    plot(xc,q_roe(1,:),'-b','LineWidth',2);
    hold all
    plot(xc,q_exact(1,:),'--k','LineWidth',2);
    ylim([0.4 1.6]);xlim([0 2]);
    legend('Lax-Friedrich method', 'Exact solution','Location','Best')
    grid on;
    title(["Time = "+num2str(Tfinal), "dx = "+num2str(dx)]) 
    ylabel('Depth')
    xlabel('x')
    hold off

    subplot(2,1,2)
    plot(xc,q_lf(2,:),'-r','LineWidth',2);
    hold all
    plot(xc,q_roe(2,:),'-b','LineWidth',2);
    hold all
    plot(xc,q_exact(2,:),'--k','LineWidth',2);
    ylim([0.1 0.4]);xlim([0 2]);
    legend('Lax-Friedrich method', 'Exact solution','Location','Best')
    grid on;
    ylabel('Discharge')
    hold off
    xlabel('x')
    saveas(fig, res_path + "/" + num2str(dx)+".png")
end

dx_iter = dx_iter + 1;
fig = figure(dx_iter);
loglog(dx_values, err_lf, 'r-'); grid on;
hold all
loglog(dx_values, err_roe, 'b-');
legend('Lax-Friedrich method', 'Godunov method','Location','northwest')
ylabel('Error')
xlabel('dx')
hold off
saveas(fig, res_path + "/" + "error.png");