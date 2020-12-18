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

% Discretization and time parameters
Tfinal = 0.5;
CFL    = 0.5;

% Domain parameters
a = 0;
b = 2;
dx = 0.001;


% Initial data set
IC = @(x) [ones(1, length(x)); (x < 1) * (-1.5)];
source = @(xc,time,u,g,h0) zeros(2, length(xc));
lf_flux = @(UL, UR, g) LF_flux(UL, UR, g);
roe_flux = @(UL, UR, g) Roe_flux(UL, UR, g);
bc  = 'open';

% Define different values of \delta x
dx_values = [0.01, 0.005, 0.001];
err_dx_lf = [];
err_dx_roe = [];

dx_iter = 0;
ref_sol_file = 'ref_sol_3.mat';
load(ref_sol_file,'ref_sol_3_dx')
load(ref_sol_file,'ref_sol_3_x')
load(ref_sol_file,'ref_sol_3')

x = {};
q = {};
for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    
    dx_iter = dx_iter + 1;
    [xc, q_lf] = FV_solver(a, b, dx, bc, IC, u, g, Tfinal, CFL, lf_flux, source);
    x{end + 1} = xc;
    q{end + 1} = q_lf;
end

fig = figure(1);
subplot(2,1,1)
plot(ref_sol_3_x,ref_sol_3(1,:),'-r','LineWidth',2, 'DisplayName', 'Reference solution');
for i = dx_iter:-1:1
    hold all
    txt = ['\Delta x = ',num2str(dx_values(i))];
    plot(x{i},q{i}(1,:),'LineWidth',2, 'DisplayName',txt);
end
hold off
ylim([0.2 1.2]);xlim([0 2]);
legend('Location','Best')
grid on;
title("Time = "+num2str(Tfinal)) 
ylabel('Depth')
xlabel('x')
hold off

subplot(2,1,2)
plot(ref_sol_3_x,ref_sol_3(2,:),'-r','LineWidth',2, 'DisplayName', 'Reference solution');
for i = dx_iter:-1:1
    hold all
    txt = ['\Delta x = ',num2str(dx_values(i))];
    plot(x{i},q{i}(2,:), 'LineWidth',2, 'DisplayName',txt);
end
hold off
xlabel('x')
ylim([-1 0.2]);xlim([0 2]);
legend('Location','Best')
grid on;
ylabel('Discharge')
xlabel('x')
hold off
saveas(fig, res_path + "/" + "LF.png")

% Roe method computation
dx_values = [0.1, 0.05, 0.001];
dx_iter = 0;
x = {};
q = {};
for dx = dx_values
    disp("Computing dx = " + num2str(dx));
    
    dx_iter = dx_iter + 1;
    [xc, q_lf] = FV_solver(a, b, dx, bc, IC, u, g, Tfinal, CFL, roe_flux, source);
    x{end + 1} = xc;
    q{end + 1} = q_lf;
end

fig = figure(2);
subplot(2,1,1)
plot(ref_sol_3_x,ref_sol_3(1,:),'--r','LineWidth',2, 'DisplayName', 'Reference solution');
for i = dx_iter:-1:1
    hold all
    txt = ['\Delta x = ',num2str(dx_values(i))];
    plot(x{i},q{i}(1,:),'LineWidth',2, 'DisplayName',txt);
end
hold off
ylim([0.2 1.2]);xlim([0 2]);
legend('Location','Best')
grid on;
title("Time = "+num2str(Tfinal)) 
ylabel('Depth')
xlabel('x')
hold off

subplot(2,1,2)
plot(ref_sol_3_x,ref_sol_3(2,:),'--r','LineWidth',2, 'DisplayName', 'Reference solution');
for i = dx_iter:-1:1
    hold all
    txt = ['\Delta x = ',num2str(dx_values(i))];
    plot(x{i},q{i}(2,:), 'LineWidth',2, 'DisplayName',txt);
end
hold off
xlabel('x')
ylim([-1 0.2]);xlim([0 2]);
legend('Location','Best')
grid on;
ylabel('Discharge')
xlabel('x')
hold off
saveas(fig, res_path + "/" + "Roe.png")