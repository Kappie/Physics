function plot_m_vs_T
  width = 0.01;
  temperatures = linspace(T_crit - width, T_crit + width, 9);
  N_values = [400];
  chi_values = [16];

  data_points = calculate_order_parameter(temperatures, chi_values, 'N_values', N_values);
  plot(temperatures, squeeze(data_points(:, 1, 1)), 'o--');

end
