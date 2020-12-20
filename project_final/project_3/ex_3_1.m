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
dx_values = [0.0008 0.001 0.0025 0.005 0.0075 0.01];
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
    
    for i=1:2
        err_lf(i, end+2-i) = norm(q_exact(i, :) - q_lf(i, :)) / length(q_lf(i, :));
        err_roe(i, end+2-i) = norm(q_exact(i, :) - q_roe(i, :)) / length(q_roe(i, :));
    end
    
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

fig = figure();
set(gcf,'position',[10,10,800,400])
hold on
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^3'];
loglog(dx_values,dx_values.^3, 'k-.', 'LineWidth',1, 'DisplayName',txt);
r = polyfit(log(dx_values), log(err_lf(1, :)), 1);
loglog(dx_values, dx_values.^r(1).*exp(r(2)), '-r','LineWidth',1, 'DisplayName', "Lax-Friedrich depth");
r = polyfit(log(dx_values), log(err_roe(1, :)), 1);
loglog(dx_values, dx_values.^r(1).*exp(r(2)), '-b','LineWidth',1, 'DisplayName', "Godunov depth");

legend('Location','northeastoutside')
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "depth_error.png");

fig = figure();
set(gcf,'position',[10,10,800,400])
hold on
txt = ['h^2'];
loglog(dx_values,dx_values.^2, 'k--', 'LineWidth',1, 'DisplayName',txt);
txt = ['h^3'];
loglog(dx_values,dx_values.^3, 'k-.', 'LineWidth',1, 'DisplayName',txt);
r = polyfit(log(dx_values), log(err_lf(2, :)), 1);
loglog(dx_values, dx_values.^r(1).*exp(r(2)), '-r', 'LineWidth',1, 'DisplayName', "Lax-Friedrich discharge");
r = polyfit(log(dx_values), log(err_roe(2, :)), 1);
loglog(dx_values, dx_values.^r(1).*exp(r(2)), '-b', 'LineWidth',1, 'DisplayName', "Godunov discharge");

legend('Location','northeastoutside')
ylabel('Error')
xlabel('dx')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');
saveas(fig, res_path + "/" + "discharge_error.png");