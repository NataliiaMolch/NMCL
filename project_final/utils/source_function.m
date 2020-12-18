function [Q] = source_function(x,t,u,g)
%SOURCE_FUNCTION Source function of the equation
h0 = @(x) 1 + 0.5*sin(pi*x);

x = x - t;

Q = [ pi/2*(u-1)*cos(pi*x);
    pi/2*cos(pi*x).*(-u + u^2+g*h0(x))];
end

