function plot_correlation_length_vs_T()
  temperature_width = 0.001;
  % 2.273435314213022
  T_pseudocrit = T_crit;
  temperatures = linspace(T_pseudocrit - temperature_width, T_pseudocrit + temperature_width, 9);

  chi_values = [9, 11, 13, 15];

  tolerances = [1e-9];

  corr_lengths = calculate_correlation_length(temperatures, chi_values, 'tolerances', tolerances);
  corr_lengths = squeeze(corr_lengths(:, :, 1));

  % returns max of each column (i.e. for each chi)
  [max_lengths, indices] = max(corr_lengths);
  temperatures(indices)

  markerplot(temperatures, corr_lengths);
  make_legend(chi_values, 'chi');
  ylabel('$\xi(\chi, t)$');
  xlabel('$T$')
  axis manual
  line([T_crit, T_crit], [0, 100000], 'LineStyle', '--');
  % line([T_pseudocrit, T_pseudocrit], [0, 100000], 'LineStyle', '--');
  % export_fig('../correlation_length_chi9-15_tol1e-9_width1e-3.pdf')

end
