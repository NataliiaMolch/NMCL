
clc
clear all
close all

% Setting : boundary condition
bc = "periodic"; % "open"
lim = 'NONE';%'NONE'; % 'MUSCLE, 'TVD'
scheme = "Roe"; %"LF"

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
dx_values = [0.05];%[0.0005 0.001 0.005 0.01]; %[ 0.01 0.05 0.1 0.5];
err_dx1 = []; err_dx2 = []; err_dx3 = []; err_dx4 = [];

dx_iter = 0;

for dx = dx_values
    disp("computing dx = " + num2str(dx))
    % Discretization
    h = dx;
    x = a : dx : b;
    xc = (a+0.5*dx):dx:(b-0.5*dx);
    N  = length(xc);
    % Initial values
    Uavg = [h0(xc); m0(xc)];
    U1 = Uavg;
%     U2 = Uavg;
%     U3 = Uavg;
%     U4 = Uavg;
    %q = (q(:,1:end-1) + q(:,2:end))/2;
    
    iter1 = 0;iter2 = 0;iter3 = 0;iter4 = 0;
    time1 = 0;
    time2 = 0; time3 = 0; time4 = 0;
    
    
    while time1 < T
        
        if sum(U1(1,:) < 0) ~= 0 
            n_complex = sum(U1(1,:) < 0);
            disp("U1 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(h) + " and iteration " + num2str(iter1))
            
        end

        k1 = CFL*h/(max(abs(U1(2,:)./U1(1,:)) + sqrt(abs(g*U1(1,:)))));
        

        % Update solution to k
%         disp(size(evalRHS(U1,g,k1,h,N,bc,lim,M, scheme)))
%         disp(size(Source_function(xc,time1,u,g,h0)))
        L = evalRHS(U1,g,k1,h,N,bc,lim,M, scheme) + Source_function(xc,time1,u,g,h0);
        Utmp = U1 + k1*L;
        
        L = evalRHS(Utmp,g,k1,h,N,bc,lim,M, scheme) + Source_function(xc,time1+k1,u,g,h0);
        Utmp = 0.75*U1 + 0.25*(Utmp + k1*L);
        
        L = evalRHS(Utmp,g,k1,h,N,bc,lim,M, scheme) + Source_function(xc,time1+k1*0.5,u,g,h0);
        U1 = U1/3 + 2/3*(Utmp + k1*L);
        
        time1 = time1 + k1;
        iter1 = iter1 + 1;
    end
    
%     while time2 < T
%         if sum(U2(1,:) < 0) ~= 0
%             n_complex = sum(U2(1,:) < 0);
%             disp("U2 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(h) + " and iteration " + num2str(iter2))
%         end
% 
%         k2 = CFL*h./(max(abs(U2(2,:)./U2(1,:)) + sqrt(abs(g*U2(1,:)))));
%         
%         RHS2 = evalRHS(U2,g,k2,h,N,bc,'MINMOD',M);
% 
%         % Update solution to k
%         U2 = U2 + RHS2;
% 
%         time2 = time2 + k2;
%         iter2 = iter2 + 1;
%     end
%     
%     while time3 < T
%         
%         if sum(U3(1,:) < 0) ~= 0
%             n_complex = sum(U3(1,:) < 0);
%             disp("U3 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(h) + " and iteration " + num2str(iter3))
%         end
% 
%         k3 = CFL*h./(max(abs(U3(2,:)./U3(1,:)) + sqrt(abs(g*U3(1,:)))));
%         
%         RHS3 = evalRHS(U3,g,k3,h,N,bc,'MUSCL',M);
% 
%         % Update solution to k
%         U3 = U3 + RHS3;
% 
%         time3 = time3 + k3;
%         iter3 = iter3 + 1;
%     end
%     
%     while time4 < T
%         
%         if sum(U4(1,:) < 0) ~= 0
%             n_complex = sum(U4(1,:) < 0);
%             disp("U4 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(h) + " and iteration " + num2str(iter4))
%         end
%         
%         k4 = CFL*h./(max(abs(U4(2,:)./U4(1,:)) + sqrt(abs(g*U4(1,:)))));
% 
%         RHS4 = evalRHS(U4,g,k4,h,N,bc,'TVB',M);
% 
%         % Update solution to k
%         U4 = U4 + RHS4;
% 
%         time4= time4 + k4;
%         iter4 = iter4 + 1;
%     end
    
%     if sum(stability_check) > 0
%         n_viol = sum(stability_check);
%         disp("Found stability violation on " + num2str(n_viol) + "iterations")
%     else
%         disp("Stable")
%     end
    
    err_dx1 = [err_dx1, norm(q_exact(xc, T) - U1,2)];
%     err_dx2 = [err_dx2, norm(q_exact(xc, T) - U2,2)];
%     err_dx3 = [err_dx3, norm(q_exact(xc, T) - U3,2)];
%     err_dx4 = [err_dx4, norm(q_exact(xc, T) - U4,2)];
%     
    dx_iter = dx_iter + 1;
    
    fig = figure(dx_iter);
    plot(xc, U1(1,:), 'r*-', xc, U1(2,:), 'b*-',x_values, q_exact_T(1,:), 'r--', x_values, q_exact_T(2,:), 'b--' ); grid on;
    legend(["Num. h", "Num. m", "Exact h", "Exact m"])
    title(["dx = "+num2str(dx)])
    saveas(fig, "ex_1/"+lim+"/"+num2str(dx)+".png")
end

dx_iter = dx_iter + 1;
fig = figure(dx_iter);
loglog(dx_values, err_dx1, 'r*-'); grid on;
title(["Error vs \delta x"])
legend([lim])
saveas(fig, "ex_1/"+lim+"/error.png")
    
