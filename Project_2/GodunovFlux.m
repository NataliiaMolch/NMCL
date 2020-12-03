function Flux = GodunovFlux(A,absA,UL,UR)

% N=size(UL, 2);
N = length(UL);
for i=1:N
    Flux(:,i) = 0.5*A(:,:,i)*(UL(:,i)+UR(:,i)) - 0.5*absA(:,:,i)*(UR(:,i)-UL(:,i));
end

end