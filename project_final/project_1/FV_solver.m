function [xc, q] = FV_solver(a, b, dx, bc, IC, u, g, Tfinal, CFL, flux, source)
% Cell centers
xc = a + 0.5*dx:dx:(b-0.5*dx);

% Cell-center values sufficient for first-order schemes
q = IC(xc);

time = 0;
iter = 0;

% Solve
while time < Tfinal    
    k = CFL*dx/(max(abs(q(2,:)./q(1,:)) + sqrt(g*q(1,:))));

    if(time + k > Tfinal)
        k = Tfinal - time;
    end

    % Applying boundary conditions to obtain extended vector
    q_ext = apply_bc_2D(q,bc);

    flux_2D = flux(q_ext(:,1:end-1),q_ext(:,2:end), g);
    
    % Compute sourse function at the current point
    if nargin < 10
        Q_ext = zeros(2, size(xc, 2));
    else
        Q_ext = source(xc,time,u,g);
    end
        
    %Scheme implementation
    q = q - k/dx*(flux_2D(:,2:end) - flux_2D(:,1:end-1)) + k*Q_ext;

    time = time + k;
    iter = iter + 1;
end

end