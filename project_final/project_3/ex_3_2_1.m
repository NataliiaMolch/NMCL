clc
clear all
close all

% Adding path to utils functions
addpath('../utils')

% Creating path for results
res_path = "2_1";
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
ICh = @(x) 1 - 0.1*sin(pi*x); 
ICm = @(x) zeros(1, length(x));
integrated_source = @(xf,time,u,g) zeros(2, length(xf)-1);
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'periodic';

% Define different values of \delta x
dx_values = [0.015, 0.01, 0.008, 0.005];
err_lf = [];
err_roe = [];
dx_iter = 0;
ref_sol_file = 'ref_sol_1.mat';
load(ref_sol_file,'ref_sol_1_dx')
load(ref_sol_file,'ref_sol_1_x')
load(ref_sol_file,'ref_sol_1')

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
        while (ref_sol_1_x(j) + 0.5 * ref_sol_1_dx) < i * dx
            j = j + 1;
        end
        ref_sol(:, end + 1) = mean(ref_sol_1(:, j_start:j), 2);
        i = i + 1;
    end
    
    err_lf(end+1) = norm(ref_sol - q_lf,2);
    err_roe(end+1) = norm(ref_sol - q_roe,2);
    
    fig = figure(dx_iter);
    subplot(2,1,1)
    plot(xc,q_lf(1,:),'-r','LineWidth',2);
    hold all
    plot(xc,q_roe(1,:),'-b','LineWidth',2);
    hold all
    plot(xc,ref_sol(1,:),'--k','LineWidth',2);
    ylim([0.8 1.2]);xlim([a b]);
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
    ylim([-0.1 0.1]);xlim([a b]);
    legend('Lax-Friedrich method', 'Godunov method','Reference solution','Location','Best')
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