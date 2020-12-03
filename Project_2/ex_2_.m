
clc
clear all
close all

% Setting : boundary condition
bc_name = "periodic"; % "open"
IC = 1; %or 2

% Domains
a = 0; b = 2.0; t0 = 0; T = 2;

% Constants
u = 0.25; CFL = 0.5; g = 1; 

% Initial conditions
switch IC
    case 1
        h0 = @(x) 1 - 0.1*sin(pi*x); 
        m0 = @(x) zeros(size(x));
    case 2
        h0 = @(x) 1 - 0.2*sin(2*pi*x); 
        m0 = @(x) 0.5*ones(size(x));
end

% Define different values of \delta x
dx_values = [0.00005 0.0001 0.0005 0.001];
dx_fine = 0.00001;
c = dx_values/dx_fine;  %compression for fine mesh
err_dx = [];

dx_iter = 0;
x_values = 0;

% Compute "exact" solution
x_fine = a : dx_fine : b;
q_fine = [h0(x_fine); m0(x_fine)];

stability_check_fine = [];
n_alpha_fine = 0;

time = 0; iter = 0;

while time < T
    if sum(q_fine(1,:) < 0) ~= 0
        n_complex = sum(q_fine(1,:) < 0);
        disp("k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx_fine) + " and iteration " + num2str(iter))
        break
    end

    k = CFL/(max(abs(q_fine(2,:)/q_fine(1,:)) + sqrt(abs(g*q_fine(1,:)))));

    if(time + k > T)
        k = T - time;
    end

    % Apply boundary conditions
    q_fine_ext = applyBC_2D(q_fine, bc_name);
    f_fine_ext = flux_func(q_fine_ext, g);

    % Compute what is necessery for LF Flux
    [Flux_fine, alpha_fine] = Lax_Friedrich_2D(q_fine_ext, f_fine_ext, g);

    % Compute sourse function at the current point
    Q_fine_ext = Source_function(x_fine,time,u,g,h0);

    % Store the violations of the stability condition
    stability_check_fine = [stability_check_fine, sum(alpha_fine > dx_fine/k)];
    n_alpha_fine = n_alpha_fine + length(alpha_fine);

    %Scheme implementation
    q_fine = q_fine - k*(Flux_fine(:,2:end) - Flux_fine(:,1:end-1)) + k*Q_fine_ext;

    time = time + k;
    iter = iter + 1;
end
    
stability_check_fine(end) = stability_check_fine(end)/n_alpha_fine;
if sum(stability_check_fine) > 0
    per_cent_viol = stability_check_fine(end);
    disp("Found stability violation in " + num2str(per_cent_viol) + " per cent of nodes")
else
    disp("Stable")
end

for dx = dx_values
    disp("computing dx = " + num2str(dx))
    % Discretization
    x = a : dx : b;
    
    % Initial values
    q = [h0(x); m0(x)];
    
    time = 0; iter = 0;
    stability_check = [];
    n_alpha = 0;
    
    while time < T

        if sum(q(1,:) < 0) ~= 0
            n_complex = sum(q(1,:) < 0);
            disp("k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx) + " and iteration " + num2str(iter))
            break
        end
        
        k = CFL/(max(abs(q(2,:)/q(1,:)) + sqrt(abs(g*q(1,:)))));
        
        if(time + k > T)
            k = T - time;
        end
                
        % Apply boundary conditions
        q_ext = applyBC_2D(q,bc_name);
        f_ext = flux_func(q_ext,g);

        % Compute what is necessery for LF Flux
        [Flux, alpha] = Lax_Friedrich_2D(q_ext, f_ext, g);

        % Compute sourse function at the current point
        Q_ext = Source_function(x,time,u,g,h0);
     
        % Store the violations of the stability condition
        stability_check = [stability_check, sum(alpha > dx/k)];
        n_alpha = n_alpha + length(alpha);

        %Scheme implementation
        q = q - k*(Flux(:,2:end) - Flux(:,1:end-1)) + k*Q_ext;
 
        time = time + k;
        iter = iter + 1;
    end
    
    
    stability_check(end) = stability_check(end)/n_alpha;
    
    if sum(stability_check) > 0
        per_cent_viol = stability_check(end);
        disp("Found stability violation in " + num2str(per_cent_viol) + " per cent of nodes")
    else
        disp("Stable")
    end
    
    dx_iter = dx_iter + 1;
    
     % Correlate mesh and fine mesh
    iter_c = 1; q_fine_compressed = zeros(size(q));
    iter = 1; compres = int32(c(dx_iter));
    while iter_c + compres <= length(q_fine)
        q_fine_compressed(:,iter) = mean(q_fine(:,iter_c:iter_c+compres), 2);
        iter_c = iter_c+compres; iter = iter + 1;
    end
    q_fine_compressed(:,end) = mean(q_fine(:,iter_c:end), 2);
    err = q - q_fine_compressed;
    err_dx = [err_dx, norm(err,2)];
    
    
    fig = figure(dx_iter);
    plot(x, q(1,:), 'r*-', x, q(2,:), 'b*-',x_fine, q_fine(1,:), 'r--', x_fine, q_fine(2,:), 'b--' ); grid on;
    legend(["Num. h", "Num. m", "Exact h", "Exact m"])
    title(["T = 2. dx = "+num2str(dx)])
    xlabel("x")
    saveas(fig, "ex_2/"+num2str(IC)+"/"+num2str(dx)+".png")
end

dx_iter = dx_iter + 1;
fig = figure(dx_iter);
loglog(dx_values, err_dx, 'k-'); grid on;
legend(["Error vs \Delta x"])
xlabel("dx")
ylabel("Error")
saveas(fig, "ex_2/"+num2str(IC)+"/"+"error"+ num2str(IC)+ ".png")
    