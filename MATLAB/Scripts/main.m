function main
  % clears cached versions of .m files (hopefully)
  rehash

  % width = 0.01;
  % T_pseudocrit = 2.274;
  % temperatures = linspace(T_pseudocrit - width, T_pseudocrit + width, 9);

  % width = 0.01;
  % temperatures = linspace(T_crit, T_crit + width, 9);

  % plot_m_vs_N_at_tcrit();
  % data_collapse_chi
  % plot_m_vs_chi_at_tcrit


  plot_correlation_length_vs_T
  % plot_m_vs_T

  % handles = plot(temperatures, squeeze(data_points(:, :, 1)));
  % set_markers(handles);
  % make_legend(chi_values, 'chi')
  % axis manual
  % line([T_crit, T_crit], [0, 100000], 'LineStyle', '--');

  % title('\chi = 64')
  % xlabel('T')
  % ylabel('\xi (correlation length)')
end
