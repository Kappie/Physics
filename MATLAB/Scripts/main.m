temperatures = [1, 2, 3];
order_parameters = ising_2d(temperatures);
save_to_file([temperatures; order_parameters]', 'testje');












  % function name = suggested_file_name()
  %   name = strcat('chi', num2str(chi), 'tolerance', num2str(tolerance), filename_suffix, '.dat')
  % end

  % data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/';
  % filename = 'chi4_adjustedreverse.dat'
