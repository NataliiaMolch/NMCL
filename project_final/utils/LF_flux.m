function Flux = LF_flux(UL, UR, g)
%LF_flux Local lax friederich flux

N=size(UL, 2);
max_eig_val1 = abs(UL(2,:)./UL(1,:))+sqrt(abs(g*UL(1,:)));
max_eig_val2 = abs(UR(2,:)./UR(1,:))+sqrt(abs(g*UR(1,:)));
alpha = max(max_eig_val1, max_eig_val2);
% alpha_LF = [alpha;alpha];

Flux = 0.5 * (flux_func(UL, g) + flux_func(UR, g)) - (0.5 * [alpha;alpha]).*(UR - UL);

end

