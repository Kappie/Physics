function plot_m_vs_chi_at_tcrit()
  temperature = T_crit;
  chi_values = [4:32];
  N_values = [100, 200, 500, 1000];

  order_parameters = squeeze(calculate_order_parameter([temperature], chi_values, 'N_values', N_values));
  markerplot(1./chi_values, order_parameters);
  make_legend(N_values, 'N')

  xlabel('$1 / \chi$');
  ylabel('$|m|$');
  % title('Order parameter at $T_c$')
  export_fig('../Plots/order_parameter_vs_1_over_chi_N100-1000.pdf')
end
