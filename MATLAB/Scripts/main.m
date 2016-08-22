function main
  clc;
  format long;
  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  T_crit = 1 / beta_crit;
  temperatures = linspace(T_crit - .1, T_crit + .1, 5);
  % temperatures = [T_crit - 0.1, T_crit - 0.1, 3]

  simulation = true;
  overwrite = false;
  plotting = false;

  min_iterations_values = [0];
  chi_values = [16];
  max_iterations_values = [60000];
  tolerance_values = [10e-30];
  chi_init_values = [2];
  tensor_initialization_values = {'adjusted'};
  traversal_order = 'reverse';

  if plotting
    figure;
    hold on;
  end

  for chi = chi_values
    for tolerance = tolerance_values
      for max_iterations = max_iterations_values
        for tensor_initialization = tensor_initialization_values
          for chi_init = chi_init_values
            for min_iterations = min_iterations_values
              file_name = suggested_file_name(chi, tolerance, max_iterations, ...
              tensor_initialization{1}, chi_init, min_iterations);

              if simulation
                order_parameters = ising_2d(temperatures, 'chi', chi, 'chi_init', chi_init, ...
                  'tolerance', tolerance, 'max_iterations', max_iterations, ...
                  'min_iterations', min_iterations, 'tensor_initialization', ...
                  tensor_initialization{1}, 'traversal_order', traversal_order);
                save_to_file([temperatures; order_parameters]', file_name, overwrite);
              end

              if plotting
                data = dlmread(file_name);
                plot(data(:, 1), data(:, 2), 'o--');
              end
            end
          end
        end
      end
    end
  end

  if plotting
    xlabel('T');
    ylabel('|m|');
    line([T_crit, T_crit], [0, 1], 'LineStyle', '--');
    hold off;
  end


  function name = suggested_file_name(chi, tolerance, max_iterations, ...
  tensor_initialization, chi_init, min_iterations)
    data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/';
    name = strcat('chi', num2str(chi), 'tolerance', num2str(tolerance, '%.0g'), 'max_iterations', ...
      num2str(max_iterations), tensor_initialization, 'chi_init', ...
      num2str(chi_init), 'min_iterations', num2str(min_iterations), '.dat');
    name = fullfile(data_dir, name);
  end

end
