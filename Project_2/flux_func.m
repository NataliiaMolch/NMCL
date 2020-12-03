function [f] = flux_func(q, g)
%FLUX_FUNC compute flux vector function of the equation
%   See the initial equation
f = [q(2,:); q(2,:).^2./q(1,:) + q(1,:).^2*g/2];
end

