function plot_m_vs_N_at_tcrit()
  temperatures = [T_crit];
  chi_values = [4, 8, 12, 16, 32, 64, 148];
  N_values = rounded_powerspace(-1/8, 0.32, 0.7, 30);

  data_points = calculate_order_parameter(temperatures, chi_values, 'N_values', N_values);
  markerplot(N_values .^ (-1/8), squeeze(data_points(1, :, :)));
  make_legend(chi_values, 'chi');

  xlabel('$N^{-1/8}$');
  ylabel('$|m|$');
  export_fig('../Plots/figure2_nishino.pdf')

end
