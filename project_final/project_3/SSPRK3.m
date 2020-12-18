function [xc, q] = SSPRK3(a, b, dx, bc, ICh, ICm, u, g, k, Tfinal, CFL, flux, integrated_source)
% Cell boundaries
xf = a:dx:b;
% Cell centers
xc = (a+0.5*dx):dx:(b-0.5*dx);

% Cell-center values sufficient for first-order schemes
N = length(xc);
q = zeros(2, N);
for j = 1:N
    q(1,j) = integral(ICh,xf(j),xf(j+1),'AbsTol',1e-14)/dx;
    q(2,j) = integral(ICm,xf(j),xf(j+1),'AbsTol',1e-14)/dx;
end

time = 0;
iter = 0;

Crec = zeros(k+1,k);
for r=-1:k-1
    Crec(r+2,:) = ReconstructWeights(k,r);
end

% Solve
while time < Tfinal        
    if sum(q(1,:) < 0) ~= 0 
        n_complex = sum(q(1,:) < 0);
        disp("U1 : k is complex " + num2str(n_complex)+ " for dx = " + num2str(dx) + " and iteration " + num2str(iter));
%         return
    end
    
    dt = CFL*dx/(max(abs(q(2,:)./q(1,:)) + sqrt(abs(g*q(1,:)))));

    if(time + dt > Tfinal)
        dt = Tfinal - time;
    end
    
    if (k == 3)
        dt = dt^(5/3);
    end
    
    %Scheme implementation
    L = evalRHS(q,bc,dx,k,Crec,flux,g) + integrated_source(xf,time,u,g);
    Utmp = q + dt*L;
    
    L = evalRHS(Utmp,bc,dx,k,Crec,flux,g) + integrated_source(xf,time+dt,u,g);
    Utmp = 0.75*q + 0.25*(Utmp + dt*L);
    
    L = evalRHS(Utmp,bc,dx,k,Crec,flux,g) + integrated_source(xf,time+0.5*dt,u,g);
    q = q/3.0 + 2.0/3.0*(Utmp + dt*L);

    time = time + dt;
    iter = iter + 1;
end

end