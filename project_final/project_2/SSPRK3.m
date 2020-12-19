function [xc, q] = MUSCL_SSPRK3(a, b, dx, bc, IC, u, g, M, Tfinal, CFL, flux, lim, source)
% Cell centers
xc = a + 0.5*dx:dx:(b-0.5*dx);
% Cell-center values sufficient for first-order schemes
q = IC(xc);

time = 0;
iter = 0;

% Solve
while time < Tfinal        
    if sum(q(1,:) < 0) ~= 0 
        n_complex = sum(q(1,:) < 0);
        disp("U1 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx) + " and iteration " + num2str(iter));
        return
    end
    
    k = CFL*dx/(max(abs(q(2,:)./q(1,:)) + sqrt(g*q(1,:))));

    if(time + k > Tfinal)
        k = Tfinal - time;
    end
        
    %Scheme implementation
    L = evalRHS(q,g,k,dx,bc,lim,M,flux) + source(xc,time,u,g);
    Utmp = q + k*L;
    
    L = evalRHS(Utmp,g,k,dx,bc,lim,M,flux) + source(xc,time+k,u,g);
    Utmp = 0.75*q + 0.25*(Utmp + k*L);
    
    L = evalRHS(Utmp,g,k,dx,bc,lim,M,flux) + source(xc,time+0.5*k,u,g);
    q = q/3 + 2/3*(Utmp + k*L);

    time = time + k;
    iter = iter + 1;
end

end