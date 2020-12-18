function RHS = evalRHS(U,bc,h,k,Crec,flux,g)

N = length(U);
U_ext = apply_bc_2D(U,bc,k);

% Obtain reconstructed states
UL = zeros(2,N+2);
UR = zeros(2,N+2);

for i = 1:N + 2
    [UL(1,i),UR(1,i)] = WENO(U_ext(1,i:(i+2*(k-1)))',k,Crec);
    [UL(2,i),UR(2,i)] = WENO(U_ext(2,i:(i+2*(k-1)))',k,Crec);
end

flux_2D = flux(UR(:, 1:end-1), UL(:, 2:end),g);

% Update solution to k
RHS =  - (flux_2D(:,2:end) - flux_2D(:,1:end-1)) / h;
return
