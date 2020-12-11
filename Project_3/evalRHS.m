% Function evaluates the RHS corresponding to a second-order MUSCL scheme
% with a mid-point rule for time-integration

function RHS = evalRHS(U,N,bc,k,h,Crec, scheme, g)

U_ext = applyBC_2D(U,bc,k);
[A, absA] = RoeSolver(U_ext,g);
% Obtain reconstructed states
UL = zeros(2,N+2);
UR = zeros(2,N+2);

for i = 1:N + 2
    [UL(1,i),UR(1,i)] = WENO(U_ext(1,i:(i+2*(k-1)))',k,Crec);
    [UL(2,i),UR(2,i)] = WENO(U_ext(2,i:(i+2*(k-1)))',k,Crec);
end

fL = flux_func(UL,g); fR = flux_func(UR,g);
% Evaluate flux
if scheme == "LF"
    flux = Lax_Friedrich_2D(U_ext,fL(:,1:end-1),fR(:,2:end),UL(:,1:end-1),UR(:,2:end),g);
elseif scheme == "Roe"
    flux = GodunovFlux(absA,fL,fR,UL,UR);
end
% Update solution to k
RHS =  - (flux(:,2:end) - flux(:,1:end-1))/h;


return
