function Flux = LF_flux(UL, UR, g)
%LF_flux Local lax friederich flux

max_eig_val1 = abs(abs(UL(2,:)/UL(1,:))+sqrt(abs(g*UL(1,:))));
max_eig_val2 = abs(abs(UR(2,:)/UR(1,:))+sqrt(abs(g*UR(1,:))));
alpha = max(max_eig_val1, max_eig_val2);
alpha_LF = [alpha;alpha];

Flux = 0.5 * (flux_func(UL, g) + flux_func(UR, g)) - alpha_LF/2.*(UR - UL);

end

