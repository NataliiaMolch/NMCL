function [F,alpha_LF] = Lax_Friedrich_2D(q,g)
%LAX_FRIEDRICH_2D Local lax friederich flux

f = flux_func(q, g);
u = q(2,:)/q(1,:);
c = sqrt(abs(g*q(1,:)));
% e_val1 = abs(u + c); e_val2 = abs(u - c);
% 
% alpha = max(max(e_val1(1:end-1), max(e_val1(2:end))), max(e_val2(1:end-1), max(e_val2(2:end))));
max_eig_val = abs(abs(u)+c);
alpha = max(max_eig_val(1:end-1), max_eig_val(2:end));
alpha_LF = [alpha;alpha];

F = (f(:,1:end-1) + f(:,2:end))/2 - alpha_LF/2.*(q(:,2:end) - q(:,1:end-1));

end

