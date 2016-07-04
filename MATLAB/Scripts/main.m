function main
  clc;
  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  T_crit = 1 / beta_crit;
  temperatures = linspace(T_crit - .1, T_crit + .1, 10);

  simulation = true;
  overwrite = true;
  plotting = true;

  min_iterations = 400;
  chi_values = [8];
  max_iterations_values = [10000];
  tolerance_values = [1];
  tensor_initialization = 'adjusted_reverse';
  % max_iterations_values = [500];

  % chi_values = [2];
  % tolerance_values = [1e-2];
  % max_iterations_values = [100];

  if simulation
    for chi = chi_values
      for tolerance = tolerance_values
        for max_iterations = max_iterations_values

          order_parameters = ising_2d(temperatures, 'chi', chi, 'chi_init', 2, ...
            'tolerance', tolerance, 'max_iterations', max_iterations, ...
            'min_iterations', min_iterations, 'tensor_initialization', tensor_initialization);

          save_to_file([temperatures; order_parameters]', ...
            suggested_file_name(chi, tolerance, max_iterations, tensor_initialization), overwrite);
        end
      end
    end
  end

  if plotting

    figure;
    hold on;
    for chi = chi_values
      for tolerance = tolerance_values
        for max_iterations = max_iterations_values
          file_name = suggested_file_name(chi, tolerance, max_iterations, tensor_initialization);
          % Want to plot columns
          data = dlmread(file_name);
          plot(data(:, 1), data(:, 2), 'o--');
        end
      end
    end

    legend(arrayfun(@num2str, chi_values, 'UniformOutput', false));


    line([T_crit, T_crit], [0, 1], 'LineStyle', '--');
    hold off;
  end

  % errors = abs(order_parameters - arrayfun(@exact_order_parameter, 1./temperatures))
  % exact_order_parameters = arrayfun(@exact_order_parameter, 1./temperatures);
  % errors = abs(order_parameters - exact_order_parameters);



  function name = suggested_file_name(chi, tolerance, max_iterations, tensor_initialization)
    data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/';
    name = strcat('chi', num2str(chi), 'tolerance', num2str(tolerance, '%.0g'), 'max_iterations', ...
      num2str(max_iterations), tensor_initialization, '.dat');
    name = fullfile(data_dir, name);
  end

end
