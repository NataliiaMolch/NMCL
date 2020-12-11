% function Flux = GodunovFlux(A,absA,UL,UR)
% 
% % N=size(UL, 2);
% N = length(UL);
% for i=1:N
%     Flux(:,i) = 0.5*A(:,:,i)*(UL(:,i)+UR(:,i)) - 0.5*absA(:,:,i)*(UR(:,i)-UL(:,i));
% end
% 
% end

function Flux = GodunovFlux(absA,fL,fR,UL,UR)


Flux = 0.5 * (fL(:,2:end) + fR(:,1:end-1));
N=size(Flux, 2);

for i=1:N
    Flux(:, i) = Flux(:, i) - 0.5*absA(:, :, i+1)*(UR(:, i)-UL(:, i+1));
end

end