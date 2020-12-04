function [F,alpha_LF] = Lax_Friedrich_2D(q,fL,fR,UL,UR,g)
%LAX_FRIEDRICH_2D Local lax friederich flux

q = q(:,2:end-1);
% u = q(2,:)/q(1,:);
% c = sqrt(abs(g*q(1,:)));
% e_val1 = abs(u + c); e_val2 = abs(u - c);
% 
% alpha = max(max(e_val1(1:end-1), max(e_val1(2:end))), max(e_val2(1:end-1), max(e_val2(2:end))));
max_eig_val = abs(abs(q(2,:)/q(1,:))+sqrt(abs(g*q(1,:))));
alpha = max(max_eig_val(1:end-1), max_eig_val(2:end));
alpha_LF = [alpha;alpha];

F = (fL + fR)/2 - alpha_LF/2.*(UR - UL);

end

