
clc
clear all
close all

% Setting : boundary condition
bc = "periodic"; % "open"
scheme = "Roe"; %"LF"
k = 2;

% Domains
a = 0; b = 2.0; t0 = 0; T = 2;

% Constants
u = 0.25; CFL = 0.5; g = 1; M = 10;

% Initial conditions
h0 = @(x) 1 + 0.5*sin(pi*x); 
m0 = @(x) u*h0(x);

%Exact solutions
q_exact = @(x,t) [h0(x - t); u*h0(x - t)];
x_values = linspace(a,b,1001);
q_exact_T = q_exact(x_values, T);

% Define different values of \delta x
dx_values = [ 0.001 0.005]; %[ 0.01 0.05 0.1 0.5];
err_dx = []; leg = [];


Crec = zeros(k+1,k);
for r=-1:k-1
    Crec(r+2,:) = ReconstructWeights(k,r);
end

for dx = dx_values
    disp("computing dx = " + num2str(dx))
    % Discretization
    h = dx;
    xf = a : dx : b;
    xc = (a+0.5*dx):dx:(b-0.5*dx);
    q_exact_arr = q_exact(xc, T);
    leg = [leg, "\Delta x = " + num2str(dx)];
    
    N  = length(xc);
    % Initial values
    U = zeros(2, N);
    for j = 1:N
        U(1,j) = integral(h0,xf(j),xf(j+1),'AbsTol',1e-14)/h;
        U(2,j) = integral(m0,xf(j),xf(j+1),'AbsTol',1e-14)/h;
    end

    iter = 0;
    time = 0;


    while time < T

        if sum(U(1,:) < 0) ~= 0 
            n_complex = sum(U(1,:) < 0);
            disp("U1 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(h) + " and iteration " + num2str(iter))

        end

        dt = CFL*h/(max(abs(U(2,:)./U(1,:)) + sqrt(abs(g*U(1,:)))));


        % Update solution to k
%         disp(size(evalRHS(U1,g,k1,h,N,bc,lim,M, scheme)))
%         disp(size(Source_function(xc,time1,u,g,h0)))
        Uold = U;

        % SSP-RK3 stage 1
        RHS = evalRHS(U,N,bc,k,h,Crec, scheme, g);
        U   = Uold + (dt/h)*RHS;

        % SSP-RK3 stage 2
        RHS = evalRHS(U,N,bc,k,h,Crec, scheme, g);
        U   = 3*Uold/4.0 + (U + (dt/h)*RHS)/4.0;

         % SSP-RK3 stage 3
        RHS = evalRHS(U,N,bc,k,h,Crec, scheme, g);
        U   = Uold/3.0 + 2.0*(U + (dt/h)*RHS)/3.0;

        time = time + dt;
        iter = iter + 1;
    end

    err_dx = [err_dx, norm(q_exact_arr - U,2)];

    fig = figure(1);

    subplot(2,1,1)
    plot(xc,U(1,:),'LineWidth',1); grid on;
    hold all

    subplot(2,1,2)
    plot(xc,U(2,:),'LineWidth',1); grid on;
    hold all

end



folder = 'ex_1/' +scheme+ '/'; 
mkdir(folder)

fig = figure(2);
loglog(dx_values, err_dx); hold all; grid on;
title(["Error vs \Delta x for different limiters"])
xlabel("\Delta x")
ylabel("Error")
saveas(fig, folder+"error_combined.png")

fig2 = figure(1);

subplot(2,1,1)
plot(x_values,q_exact_T(1,:),'--k','LineWidth',2);
xlim([0 2]);
title("Depth")
legend(leg,'Location','Best')
grid on;
hold off

subplot(2,1,2)
plot(x_values,q_exact_T(2,:),'--k','LineWidth',2);
xlim([0 2]);
title("Discharge")
legend(leg,'Location','Best')
grid on;
hold off
dx = dx_values(1);
filename = folder+num2str(dx)+".png";
saveas(fig, filename)