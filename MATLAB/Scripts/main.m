function main
  clc;
  format long;
  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  T_crit = 1 / beta_crit;

  % max_iterations = 1e6;
  % tolerance = 1e-7;
  temperature = T_crit + 0.01;
  chi = 8;
  N = 250;

  result = ising_2d([temperature], 'chi', chi, 'N', N)

  % Experiment 1
  % Plot magnetization vs temperature around critical point.

  % temperatures = linspace(T_crit - .01, T_crit + .01, 9);
  % chi_values = [2, 4, 8, 16, 32, 48, 64, 80];
  % legend_labels = arrayfun(@(chi) ['chi = ' num2str(chi)], chi_values, 'UniformOutput', false);
  %
  % dataset = zeros(numel(temperatures), numel(chi_values));
  %
  % for i = 1:numel(chi_values)
  %   order_parameters = ising_2d(temperatures, 'chi', chi_values(i), ...
  %     'tolerance', tolerance, 'max_iterations', max_iterations);
  %   dataset(:,i) = order_parameters;
  % end
  %
  % figure
  % plot(temperatures, dataset, '--o')
  % legend(legend_labels, 'Location', 'southwest')
  % line([T_crit, T_crit], [0, 1], 'LineStyle', '--');
  % xlabel('T')
  % ylabel('|m|')
  % title('Magnetization of 2D Ising model at tolerance = 1e-7')

  % Experiment 2
  % Plot magnetization versus 1/chi for high chi_values.

  % chi_values = [2, 4, 8, 16, 32, 64, 80, 96, 112];
  % tolerance_values = [1e-4, 1e-5, 1e-6, 1e-7];
  % temperature = T_crit
  % legend_labels = arrayfun(@(tolerance) ['tolerance = ' num2str(tolerance, '%.0e')], tolerance_values, 'UniformOutput', false);
  %
  % magnetizations = zeros(numel(chi_values), numel(tolerance_values));
  %
  % for t = 1:numel(tolerance_values)
  %   for c = 1:numel(chi_values)
  %     profile on
  %     result = ising_2d([temperature], 'chi', chi_values(c), 'tolerance', ...
  %       tolerance_values(t), 'max_iterations', max_iterations);
  %     profsave
  %     magnetizations(c, t) = result(1);
  %   end
  % end

  % figure
  % plot(1./chi_values, magnetizations, 'o--')
  % legend(legend_labels);
  % xlabel('1 / chi')
  % ylabel('|m|')
  % title('Magnetization of 2D Ising model at T_{crit}')

end
