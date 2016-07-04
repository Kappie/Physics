function plot_order_parameter(temperatures, order_parameter_set)
beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
T_crit = 1 / beta_crit;

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot

hold on
for i = 1:length(order_parameter_set)
  plot(temperatures, order_parameter_set{i}, 'Marker', 'o', 'LineStyle', 'none');
end

% Create xlabel
xlabel({'T'});

% Create title
title({'Order parameter of 2D Ising model'});

% Create ylabel
ylabel({'m'});

box(axes1,'on');
% Create arrow
line([T_crit, T_crit], [0, 1], 'LineStyle', '--');

end
