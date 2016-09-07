function plot_m_vs_N_at_tcrit()
  temperature = T_crit;
  chi_values = [4, 10, 148];
  N_values = fliplr( arrayfun(@round, linspace(0.3, 0.7, 30) .^ -8) );
  MARKERS = markers();

  legend_labels = arrayfun(@(chi) ['chi = ' num2str(chi)], chi_values, 'UniformOutput', false);

  dataset = zeros(numel(N_values), numel(chi_values));

  % figure
  % hold on

  for c = 1:numel(chi_values)
    for n = 1:numel(N_values)
      order_parameters = ising_2d([temperature], 'chi', chi_values(c), 'N', N_values(n));
      dataset(n,c) = order_parameters(1);
    end
    save('nishino_fig2.mat', 'chi_values', 'N_values', 'dataset')
    % plot(N_values, dataset(:, c), MARKERS(c))
  end

  % legend(legend_labels)
  % xlabel('N')
  % ylabel('|m|')
  % title('Order parameter at T_c')
  % hold off
end
