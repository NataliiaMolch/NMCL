function [F,alpha_LF] = Lax_Friedrich_2D(q,fL,fR,UL,UR,g)
%LAX_FRIEDRICH_2D Local lax friederich flux

q = q(:,2:end-1);

max_eig_val = abs(abs(q(2,:)/q(1,:))+sqrt(abs(g*q(1,:))));
alpha = max(max_eig_val(1:end-1), max_eig_val(2:end));
alpha_LF = [alpha;alpha];

% max_eig_val1 = abs(abs(UL(2,:)/UL(1,:))+sqrt(abs(g*UL(1,:))));
% max_eig_val2 = abs(abs(UR(2,:)/UR(1,:))+sqrt(abs(g*UR(1,:))));
% alpha = max(max_eig_val1, max_eig_val2);
% alpha_LF = [alpha;alpha];

F = (fL + fR)/2 - alpha_LF/2.*(UR + UL);

end

