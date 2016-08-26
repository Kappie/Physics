function main
  clc;
  format long;
  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  T_crit = 1 / beta_crit;

  max_iterations = 1e6;
  tolerance = 1e-7;

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

  chi_values = [2, 4, 8, 16, 32, 64, 80, 96, 112];
  tolerance_values = [1e-4, 1e-5, 1e-6, 1e-7];
  temperature = T_crit 
  legend_labels = arrayfun(@(tolerance) ['tolerance = ' num2str(tolerance, '%.0e')], tolerance_values, 'UniformOutput', false);

  magnetizations = zeros(numel(chi_values), numel(tolerance_values));

  for t = 1:numel(tolerance_values)
    for c = 1:numel(chi_values)
      profile on
      result = ising_2d([temperature], 'chi', chi_values(c), 'tolerance', ...
        tolerance_values(t), 'max_iterations', max_iterations);
      profsave
      magnetizations(c, t) = result(1);
    end
  end

  % figure
  % plot(1./chi_values, magnetizations, 'o--')
  % legend(legend_labels);
  % xlabel('1 / chi')
  % ylabel('|m|')
  % title('Magnetization of 2D Ising model at T_{crit}')





  % chi_values = [2, 4, 8, 16];
  % max_iterations_values = [1e5];
  % tolerance_values = [1e-8];
  % min_iterations_values = [0];
  % chi_init_values = [2];
  % tensor_initialization_values = {'spin-up'};
  % traversal_order = 'standard';

  % simulation = true;
  % overwrite = true;
  % plotting = true;

  % if plotting
  %   figure;
  %   hold on;
  % end

  % for chi = chi_values
  %   for tolerance = tolerance_values
  %     for max_iterations = max_iterations_values
  %       for tensor_initialization = tensor_initialization_values
  %         for chi_init = chi_init_values
  %           for min_iterations = min_iterations_values
  %             file_name = suggested_file_name(chi, tolerance, max_iterations, ...
  %             tensor_initialization{1}, chi_init, min_iterations);
  %
  %             if simulation
  %               order_parameters = ising_2d(temperatures, 'chi', chi, 'chi_init', chi_init, ...
  %                 'tolerance', tolerance, 'max_iterations', max_iterations, ...
  %                 'min_iterations', min_iterations, 'tensor_initialization', ...
  %                 tensor_initialization{1}, 'traversal_order', traversal_order);
  %               save_to_file([temperatures; order_parameters]', file_name, overwrite);
  %             end
  %
  %             if plotting
  %               data = dlmread(file_name);
  %               plot(data(:, 1), data(:, 2), 'o--', 'DisplayName', ['chi = ' num2str(chi)]);
  %             end
  %           end
  %         end
  %       end
  %     end
  %   end
  % end
  %
  % if plotting
  %   legend('-DynamicLegend');
  %   xlabel('T');
  %   ylabel('|m|');
  %   line([T_crit, T_crit], [0, 1], 'LineStyle', '--');
  %   hold off;
  % end


  function name = suggested_file_name(chi, tolerance, max_iterations, ...
  tensor_initialization, chi_init, min_iterations)
    data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/';
    name = strcat('chi', num2str(chi), 'tolerance', num2str(tolerance, '%.0g'), 'max_iterations', ...
      num2str(max_iterations), tensor_initialization, 'chi_init', ...
      num2str(chi_init), 'min_iterations', num2str(min_iterations), '.dat');
    name = fullfile(data_dir, name);
  end

end
