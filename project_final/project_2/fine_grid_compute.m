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
bc  = 'open';
limiter = "TVB";

% Define different values of \delta x
dx =  0.00025;

[xc, ref_sol_3] = SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, lf_flux, limiter, source);

ref_sol_3_dx = dx;
ref_sol_3_x = xc;

save("ref_sol_3.mat", 'ref_sol_3', 'ref_sol_3_dx', 'ref_sol_3_x')