function U_bc = applyBC_2D(Ui,bc_name,m)
%APPLYBC takes array of values of the function in the nodes and adds additional
%values at the beggining and in the end that correspond either to periodic
%or open boundary conditions
% bc_name = either "periodic" or "open"
% works with 2D arrays
switch bc_name
    case "periodic"
        U_bc = [Ui(:,end-m+1:end) , Ui, Ui(:,1:m)];
    case "open"
        U_bc = [repmat(Ui(:,1),1,m), Ui, repmat(Ui(:,end),1,m)];
end

