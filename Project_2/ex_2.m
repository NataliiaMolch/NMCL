
clc
clear all
close all

% Setting : boundary condition
bc = "periodic"; % "open"
scheme = "LF";
% IC = 3; %or 2
limiters = ["NONE" "MUSCL" "MINMOD" "TVB"];

% Domains
a = 0; b = 2.0; t0 = 0; T = 2;

% Constants
u = 0.25; CFL = 0.5; g = 1; M = 10;

% Initial conditions

for IC = [1 2 3]
    switch IC
        case 1
            h0 = @(x) 1 - 0.1*sin(pi*x); 
            m0 = @(x) zeros(size(x));
        case 2
            h0 = @(x) 1 - 0.2*sin(2*pi*x); 
            m0 = @(x) 0.5*ones(size(x));
        case 3
            h0 = @(x) ones(size(x));
            m0 = @(x) -1.5*(x <= 1) + 0*(x>1);

    end

    % Define different values of \delta x
    dx_values = [0.0005 0.001 0.005 0.01];
    dx_fine = 0.0001;
    c = dx_values/dx_fine;  %compression for fine mesh
    q_exact_T = 0;

    folder = "ex_2/"+scheme+"/"+num2str(IC)+"/";
    mkdir(folder);

    % Compute "exact" solution
    x_fine = (a+0.5*dx_fine):dx_fine:(b-0.5*dx_fine);
    U_fine = [h0(x_fine); m0(x_fine)];
    N  = length(x_fine);

    time = 0; iter = 0;
    lim = "TVB";
    err_all = [];

    while time < T
        if sum(U_fine(1,:) < 0) ~= 0
            n_complex = sum(U_fine(1,:) < 0);
            disp("k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx_fine) + " and iteration " + num2str(iter))
            break
        end

        k = CFL*dx_fine/(max(abs(U_fine(2,:)./U_fine(1,:)) + sqrt(abs(g*U_fine(1,:)))));

        if(time + k > T)
            k = T - time;
        end

        % Apply boundary conditions
        h = dx_fine;
        L = evalRHS(U_fine,g,k,h,N,bc,lim,M, scheme);
        Utmp = U_fine + k*L;

        L = evalRHS(Utmp,g,k,h,N,bc,lim,M, scheme);
        Utmp = 0.75*U_fine + 0.25*(Utmp + k*L);

        L = evalRHS(Utmp,g,k,h,N,bc,lim,M, scheme);
        U_fine = U_fine/3 + 2/3*(Utmp + k*L);

        time = time + k;
        iter = iter + 1;
    end

    q_exact_T = U_fine;

    for lim = limiters
        err_dx = []; 
        dx_iter = 0;

        time = 0; iter = 0;
        for dx = dx_values
            disp("computing dx = " + num2str(dx))
            % Discretization
            x = a : dx : b;
            xc = (a+0.5*dx):dx:(b-0.5*dx);
            N  = length(xc);
            % Initial values
            U1 = [h0(xc); m0(xc)];

            time = 0; iter = 0;
            n_alpha = 0;

            while time < T

                if sum(U1(1,:) < 0) ~= 0
                    n_complex = sum(U1(1,:) < 0);
                    disp("k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx) + " and iteration " + num2str(iter))
                    break
                end

                k = CFL*dx/(max(abs(U1(2,:)./U1(1,:)) + sqrt(abs(g*U1(1,:)))));

                if(time + k > T)
                    k = T - time;
                end

                h = dx;
                L = evalRHS(U1,g,k,h,N,bc,lim,M, scheme);
                Utmp = U1 + k*L;

                L = evalRHS(Utmp,g,k,h,N,bc,lim,M, scheme);
                Utmp = 0.75*U1 + 0.25*(Utmp + k*L);

                L = evalRHS(Utmp,g,k,h,N,bc,lim,M, scheme);
                U1 = U1/3 + 2/3*(Utmp + k*L);

                time = time + k;
                iter = iter + 1;
            end

            dx_iter = dx_iter + 1;

             % Correlate mesh and fine mesh
            iter_c = 1; q_fine_compressed = zeros(size(U1));
            iter = 1; compres = int32(c(dx_iter));
            while iter_c + compres <= length(U_fine)
                q_fine_compressed(:,iter) = mean(U_fine(:,iter_c:iter_c+compres), 2);
                iter_c = iter_c+compres; iter = iter + 1;
            end
            q_fine_compressed(:,end) = mean(U_fine(:,iter_c:end), 2);
            err = U1 - q_fine_compressed;
            err_dx = [err_dx, norm(err,2)];
            
            if dx == dx_values(1)
                fig = figure(1);

                subplot(2,1,1)
                plot(xc,U1(1,:),'LineWidth',1); grid on;
                hold all

                subplot(2,1,2)
                plot(xc,U1(2,:),'LineWidth',1); grid on;
                hold all
                
            end

        end

        err_all = [err_all; err_dx];


    end

    fig2 = figure(2);

    loglog(dx_values, err_all(1,:), 'r*-',dx_values, err_all(2,:), 'm*-', ...
        dx_values, err_all(3,:), 'g*-',dx_values, err_all(4,:), 'b*-' )
    legend(["NONE" "MUSCL" "MINMOD" "TVB"])
    title(["Error vs \Delta x for different limiters"])
    xlabel("\Delta x, \Delta x = "+num2str(dx(1)))
    ylabel("Error")
    saveas(fig2, folder+"error_combined.png")

    fig1  = figure(1);
    leg = [limiters, "Exact"];

    subplot(2,1,1)
    plot(x_fine,q_exact_T(1,:),'--k','LineWidth',2);
    xlim([0 2]);
    legend(leg,'Location','Best')
    grid on;
    title("Depth, \Delta x = "+num2str(dx_values(1)))
    hold off

    subplot(2,1,2)
    plot(x_fine,q_exact_T(2,:),'--k','LineWidth',2);
    xlim([0 2]);
    title("Discharge, \Delta x = "+num2str(dx_values(1)))
    legend(leg,'Location','Best')
    grid on;
    hold off
    dx = dx_values(1);
    filename = folder+num2str(dx)+".png";
    saveas(fig1, filename)
end