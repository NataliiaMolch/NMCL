function res = plot_support(a, b, dx, xc, Tfinal, figure_num, x_exact, q_exact, limiters, q, flux_name, res_path, y1_min, y1_max, y2_min, y2_max)
fig = figure(figure_num);
subplot(2,1,1)
plot(x_exact, q_exact(1, :),'--k','LineWidth',1, 'DisplayName', 'Exact');
for i = 1:length(limiters)
    hold all
    txt = num2str(limiters(i));
    plot(xc,q{i}(1, :), 'LineWidth',1, 'DisplayName',txt);
end
hold off
ylim([y1_min y1_max]);xlim([a b]);
legend('Location','northeastoutside')
grid on;
title([flux_name + " flux","Time = "+num2str(Tfinal), "\Delta x = "+num2str(dx)]) 
ylabel('Depth')
xlabel('x')
hold off

subplot(2,1,2)
plot(x_exact, q_exact(2, :),'--k','LineWidth',1, 'DisplayName', 'Exact');
for i = 1:length(limiters)
    hold all
    txt = num2str(limiters(i));
    plot(xc,q{i}(2, :), 'LineWidth',1, 'DisplayName',txt);
end
hold off
ylim([y2_min y2_max]);xlim([a b]);
legend('Location','northeastoutside')
grid on;
ylabel('Discharge')
xlabel('x')
hold off
saveas(fig, res_path + "/" + flux_name + "_" + num2str(dx)+".png")
end