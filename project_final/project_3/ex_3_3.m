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
k = 2;

% Discretization and time parameters
Tfinal = 2;
CFL    = 0.5;

% Domain parameters
a = 0;
b = 2;

% Initial data set
ICh = @(x) ones(1, length(x)); 
ICm = @(x) (x < 1) * (-1.5);
integrated_source = @(xf,time,u,g) zeros(2, length(xf)-1);
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'periodic';

% Define different values of \delta x
dx_values = [0.01, 0.005, 0.001];
err_lf = [];
err_roe = [];
dx_iter = 0;
ref_sol_file = 'ref_sol_3.mat';
load(ref_sol_file,'ref_sol_3_dx')
load(ref_sol_file,'ref_sol_3_x')
load(ref_sol_file,'ref_sol_3')

for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    dx_iter = dx_iter + 1;
    
    xc = a+0.5*dx:dx:b-0.5*dx;

    [xc, q_lf] = SSPRK3(a, b, dx, bc, ICh, ICm, u, g, k, Tfinal, CFL, lf_flux, integrated_source);
    [xc, q_roe] = SSPRK3(a, b, dx, bc, ICh, ICm, u, g, k, Tfinal, CFL, roe_flux, integrated_source);
    
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
    
    for i=1:2
        err_lf(i, end+2-i) = norm(q_lf(i, :) -  ref_sol(i, :)) / length(q_lf(i, :));
        err_roe(i, end+2-i) = norm(q_roe(i, :) - ref_sol(i, :)) / length(q_roe(i, :));
    end
    
    fig = figure(dx_iter);
    subplot(2,1,1)
    plot(xc,q_lf(1,:),'-r','LineWidth',2);
    hold all
    plot(xc,q_roe(1,:),'-b','LineWidth',2);
    hold all
    plot(xc,ref_sol(1,:),'--k','LineWidth',2);
    ylim([0.2 1.2]);xlim([a b]);
    legend('Lax-Friedrich method', 'Godunov method','Reference solution','Location','Best')
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
    plot(xc,ref_sol(2,:),'--k','LineWidth',2);
    ylim([-1 0.2]);xlim([a b]);
    legend('Lax-Friedrich method', 'Godunov method','Reference solution','Location','Best')
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