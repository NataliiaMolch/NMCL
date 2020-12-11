function [um,up] = WENO(uloc,m,Crec)
% Purpose: Compute the left and right cell interface values using an WENO
%          approach based on 2m-1 long vectors uloc with cell
% Set WENO parameters

vareps = 1e-6;
upl = zeros(m,1); uml = zeros(m,1);

% Treat special case of m=1 - no stencil to select
if (m==1)
    um = uloc(1); up = uloc(1);
    return
elseif (m==2)
    assert(length(uloc) == 3);
    betar(1) = (uloc(3) - uloc(2))^2;
    betar(2) = (uloc(2) - uloc(1))^2;
    dw = [2.0/3.0,1.0/3.0];
    
elseif (m==3)
    assert(length(uloc) == 5);
    betar(1) = 13.0/12.0*(uloc(3) - 2*uloc(4) + uloc(5))^2 +...
               1.0/4.0*(3*uloc(3) - 4*uloc(4) + uloc(5))^2;
    betar(2) = 13.0/12.0*(uloc(2) - 2*uloc(3) + uloc(4))^2 +...
               1.0/4.0*(uloc(2) - uloc(4))^2;
    betar(3) = 13.0/12.0*(uloc(1) - 2*uloc(2) + uloc(3))^2 +...
               1.0/4.0*(uloc(1) - 4*uloc(2) + 3*uloc(3))^2;
    dw = [3.0/10.0,3.0/5.0,1.0/10.0];
    
    
end

% Compute um and up based on different stencils 
for r=0:m-1
    umh = uloc(m-r+[0:m-1]);
    upl(r+1) = Crec(r+2,:)*umh; uml(r+1) = Crec(r+1,:)*umh;
end


% Compute alpha weights - classic WENO
alphap = dw./(vareps+betar).^2;
alpham = flipud(dw)./(vareps+betar).^2;


% Compute nonlinear weights and cell interface values
um=alpham*uml/sum(alpham); up=alphap*upl/sum(alphap);

return