clf
clear all
close all

scheme = "LF";
file = 'ex_1/' +scheme+ '/vars.mat';
folder = 'ex_1/' +scheme+ '/';
load(file);

dx_values = [0.0005 0.001 0.005 0.01];

fig = figure(1);
loglog(dx_values, err_all(1,:), 'r*-',dx_values, err_all(2,:), 'm*-', ...
    dx_values, err_all(3,:), 'g*-',dx_values, err_all(4,:), 'b*-' )
legend(["NONE" "MUSCL" "MINMOD" "TVB"])
title(["Error vs \Delta x for different limiters"])
xlabel("\Delta x")
ylabel("Error")
saveas(fig, folder+"error_combuned.png")

fig = figure(2);

leg = ["NONE" "MUSCL" "MINMOD" "TVB", "Exact"];
subplot(2,1,1)
x = 0:dx_values(1):2;
plot(x,q_all(1,1,:),'-r','LineWidth',1);
x = 0:dx_values(2):2;
plot(x,q_all(2,1,:),'-m','LineWidth',1);
x = 0:dx_values(3):2;
plot(x,q_all(3,1,:),'-g','LineWidth',1);
x = 0:dx_values(4):2;
plot(x,q_all(4,1,:),'-b','LineWidth',1);
hold all
plot(xc,q_exact_arr(1,:),'--k','LineWidth',2);
xlim([0 2]);
title("Depth")
legend(leg,'Location','Best')
grid on;
title([scheme + " " + lim + " : Time = "+num2str(time1), "dx = "+num2str(dx)])
hold off

subplot(2,1,2)
x = 0:dx_values(1):2;
plot(x,q_all(1,2,:),'-r','LineWidth',1);
x = 0:dx_values(2):2;
plot(x,q_all(2,2,:),'-m','LineWidth',1);
x = 0:dx_values(3):2;
plot(x,q_all(3,2,:),'-g','LineWidth',1);
x = 0:dx_values(4):2;
plot(x,q_all(4,2,:),'-b','LineWidth',1);
hold all
plot(xc,q_exact_arr(2,:),'--k','LineWidth',2);
xlim([0 2]);
title("Discharge")
legend(leg,'Location','Best')
grid on;
hold off
filename = folder+num2str(dx)+".png";
saveas(fig, filename)