beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
T_crit = 1 / beta_crit;

temperatures = linspace(T_crit - 0.1, T_crit + 0.1, 20)
order_parameters = ising_2d(temperatures)
% save_to_file([temperatures; order_parameters]', 'testje');

plot_order_parameter(temperatures, order_parameters);












  % function name = suggested_file_name()
  %   name = strcat('chi', num2str(chi), 'tolerance', num2str(tolerance), filename_suffix, '.dat')
  % end

  % data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/';
  % filename = 'chi4_adjustedreverse.dat'
