function [Q] = Source_function(x,t,u,g,h0)
%SOURCE_FUNCTION Source function of the equation
% u = 0.25;
% h0 = @(x) 1 + 0.5*sin(pi*x);

x = x - t;

Q = [ pi/2*(u-1)*cos(pi*x);
    pi/2*cos(pi*x).*(-u + u^2+g*h0(x))];
end

