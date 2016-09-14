function main
  % width = 0.01;
  % T_pseudocrit = 2.274;
  % temperatures = linspace(T_pseudocrit - width, T_pseudocrit + width, 9);

  width = 0.01;
  temperatures = linspace(T_crit, T_crit + width, 9);
  % temperatures = [T_crit];

  chi_values = [8, 16];
  N_values = [800, 3000];

  data_points = calculate_correlation_length(temperatures, chi_values, N_values);

  handles = plot(temperatures, squeeze(data_points(:, :, 1)));
  set_markers(handles);
  make_legend(chi_values, 'chi')
  axis manual
  line([T_crit, T_crit], [0, 100000], 'LineStyle', '--');

  % title('\chi = 64')
  % xlabel('T')
  % ylabel('\xi (correlation length)')
end
