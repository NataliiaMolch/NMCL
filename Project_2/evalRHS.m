% Function evaluates the RHS corresponding to a second-order MUSCL scheme
% with a mid-point rule for time-integration

function RHS = evalRHS(U,g,k,h,N,bc,limiter,M, scheme)

% Need to extend by 2 ghost cells on either side
% Evaluate flux
if scheme == "LF"
    U_ext = applyBC_2D(U,bc,1);
    flux = Lax_Friedrich_2D(U_ext,g);
elseif scheme == "Roe"
    U_ext = applyBC_2D(U,bc,2);
    [A, absA] = RoeSolver(U_ext,g);

    % Obtain limited slope for N+2 cells
    dU      = zeros(2,N+2);
    dU(1,:) = SlopeLimiter(U_ext(1,:),limiter, M, h);
    dU(2,:) = SlopeLimiter(U_ext(2,:),limiter, M, h);

    % Obtain cell solution at k/2 with f'(u) = A
    Unph = [U_ext(:,2) - 0.5*k/h*A(:,:,1)*dU(:,1)];
    for i = 2:(length(U_ext)-2)
        Unph = [Unph,U_ext(:,i+1) - 0.5*k/h*A(:,:,i)*dU(:,i)];
    end

    % Obtain interface values at k/2
    UL = Unph(:,1:end-1) + 0.5*dU(:,1:end-1);
    UR = Unph(:,2:end)   - 0.5*dU(:,2:end);
    
    flux = GodunovFlux(A,absA,UL,UR);
end
% Update solution to k
RHS =  - (flux(:,2:end) - flux(:,1:end-1))/h;


return
