
clc
clear all
close all

% Setting : boundary condition
bc = "periodic"; % "open"
% lim = 'NONE';%'NONE'; % 'MUSCLE, 'TVD'
limiters = ["NONE" "MUSCLE" "MINMOD" "TVD"];
scheme = "LF"; %"LF"

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
dx_values = [0.0005 0.001 0.005 0.01]; %[ 0.01 0.05 0.1 0.5];


for lim = limiters
    err_dx1 = []; err_all = []; q_all = [];
    dx_iter = 0;
    iter_lim = 0;
    
    for dx = dx_values
        disp("computing dx = " + num2str(dx))
        % Discretization
        h = dx;
        x = a : dx : b;
        xc = (a+0.5*dx):dx:(b-0.5*dx);
        q_exact_arr = q_exact(xc, T);

        N  = length(xc);
        % Initial values
        Uavg = [h0(xc); m0(xc)];
        U1 = Uavg;

        iter1 = 0;
        time1 = 0;

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
        
        err_dx1 = [err_dx1, norm(q_exact_arr - U1,2)];
        
        if dx == min(dx_values)
            err_all = [err_all, err_dx];
            q_all = [q_all, U1];
        end
        dx_iter = dx_iter + 1;


        fig = figure(dx_iter);
        subplot(2,1,1)
        plot(xc,U1(1,:),'-r','LineWidth',2);
        hold all
        plot(xc,q_exact(1,:),'--k','LineWidth',2);
        ylim([0.8 1.2]);xlim([0 2]);
        legend('Numerical Depth','Exact Depth','Location','Best')
        grid on;
        title(["Time = "+num2str(time), "dx = "+num2str(dx)])
        hold off

        subplot(2,1,2)
        plot(xc,U1(2,:),'-r','LineWidth',2);
        hold all
        plot(xc,q_exact(2,:),'--k','LineWidth',2);
        ylim([-0.2 0.1]);xlim([0 2]);
        legend('Numerical Discharge','Exact Discharge','Location','Best')
        grid on;
        hold off
        folder = "ex_1/"+scheme+lim+"/";
        mkdir(folder)
        saveas(fig, folder+num2str(dx)+".png")
    end
    
    dx_iter = dx_iter + 1;
    fig = figure(dx_iter);
    loglog(dx_values, err_dx1, 'r*-'); grid on;
    title(["Error vs \Delta x"])
    xlabel("\Delta x")
    ylabel("Error")
    legend([lim])
    saveas(fig, folder+"/error.png")
    
    iter_lim  = iter_lim +1;
end
    
save("ex_1/"+scheme+lim+"/vars.mat", 'err_all', 'q_all', 'q_exact_arr')

