function fit_power_law
  load('nishino_fig2.mat', 'N_values', 'chi_values', 'dataset');
  N_values(numel(N_values))
  % slope = logfit(N_values, dataset(:, 3), 'loglog', 'skipBegin', 19)

  % handles = plot(N_values .^ -0.125, dataset);
  % set_markers(handles);

  % legend({'\chi = 4', '\chi = 10', '\chi = 148'}, 'location', 'SouthEast')
  % xlabel('N^{-1/8}')
  % ylabel('|m|')
  % title('Order parameter of Ising model T = T_c.')
  % export_fig '../Plots/figure2_nishino.pdf'
end
