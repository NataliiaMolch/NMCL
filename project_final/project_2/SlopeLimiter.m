% function dU = SlopeLimiter(U,limiter, M,h)
% 
% dUL = U(2:end-1)-U(1:end-2);
% dUR = U(3:end)-U(2:end-1);
% switch limiter
%     case 'NONE'
%         dU = 0*dUL;
%     case 'MINMOD'
%         dU = minmod([dUL',dUR'])';
%     case 'MUSCL'
%         dU = minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
%     case 'TVB'
%         a1 = 0.5*(dUL'+dUR');
%         dU = a1.*(abs(a1) <= M*h^2) + (abs(a1) > M*h^2).*minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
%         disp(size(dU))
%     otherwise
%         error('Unknown limiter function requested!')
%         
% end
%         
% 
% return
function dU = SlopeLimiter(U,limiter, M,h)

dUL = U(2:end-1)-U(1:end-2);
dUR = U(3:end)-U(2:end-1);
switch limiter
    case 'NONE'
        dU = 0*dUL;
    case 'MINMOD'
        dU = minmod([dUL',dUR'])';
    case 'MUSCL'
        dU = minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
    case 'TVB'
        a1 = dUL';
        dU = a1'.*(abs(a1) <= M*h^2)' + (abs(a1) > M*h^2)'.*minmod([dUL', dUR'])';
%         a1 = 0.5*(dUL'+dUR');
%         dU = a1'.*(abs(a1) <= M*h^2)' + (abs(a1) > M*h^2)'.*minmod([0.5*(dUL'+dUR'),2*dUL',2*dUR'])';
%         disp(limiter)
%         disp(size(dU))
    otherwise
        error('Unknown limiter function requested!')
        
end
        

return