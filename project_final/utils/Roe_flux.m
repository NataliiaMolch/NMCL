function Flux = Roe_flux(UL,UR,g)

N=size(UL, 2);
% Flux = 0.5 * (flux_func(UL, g) + flux_func(UR, g));
[absA, A] = Roe_matrix(UL, UR, g);

for i=1:N
    Flux(:, i) = 0.5 * A(:,:,i) *(UR(:, i)+UL(:, i))  - 0.5*absA(:, :, i)*(UR(:, i)-UL(:, i));
end

end