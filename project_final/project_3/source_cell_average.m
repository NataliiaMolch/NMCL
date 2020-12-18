function source = source_cell_average(xf,time,u,g)
dx = xf(2) - xf(1);
xf = pi*(xf - time);
source = [(u-1)/2*(sin(xf(2:end)) - sin(xf(1:end-1)));
          0.5*(sin(xf(2:end)) - sin(xf(1:end-1)))*(u^2-u+g)+0.125*g*(sin(xf(2:end)).^2 - sin(xf(1:end-1)).^2)] / dx;
end



