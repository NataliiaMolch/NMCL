function [absA, A] = Roe_matrix(ql, qr, g)

ql = ql(1,:).^(-0.5) .* ql;
qr = qr(1,:).^(-0.5) .* qr;
mean_z = 0.5 * (ql + qr);
mean_z1_2 = 0.5 * (ql(1,:) .* ql(1,:) + qr(1,:) .* qr(1,:));
N = length(ql(1, :));
A = reshape([zeros(1, N); 
    mean_z1_2 * g - mean_z(2, :).*mean_z(2, :) ./ (mean_z(1, :).*mean_z(1, :)); 
    ones(1, N); 
    2*mean_z(2, :) ./ mean_z(1, :)], 2, 2, []);

a = mean_z(2, :) ./ mean_z(1, :);
b = sqrt(mean_z1_2 *g);
c = a./b;
divided_sum_ab = (c + 1);
divided_diff_ab = (c - 1);

abs_lambda1 = abs(a+b);
abs_lambda2 = abs(a-b);

absA = reshape(0.5 * [abs_lambda2.*divided_sum_ab - abs_lambda1.*divided_diff_ab; 
    (abs_lambda2 - abs_lambda1) .* (divided_sum_ab .* (b - a));
    abs_lambda1-abs_lambda2;
    abs_lambda1 .* divided_sum_ab - abs_lambda2 .* divided_diff_ab], 2, 2, []);
end