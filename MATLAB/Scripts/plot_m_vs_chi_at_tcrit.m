function plot_m_vs_chi_at_tcrit()
  temperature = T_crit;
  chi_values = [4:4:64];
  N_values = [100, 1000, 10000];

  legend_labels = arrayfun(@(N) ['N = ' num2str(N)], N_values, 'UniformOutput', false);
  MARKERS = markers;

  dataset = zeros(numel(chi_values), numel(N_values));

  figure
  hold on

  for n = 1:numel(N_values)
    for c = 1:numel(chi_values)
      pause(0.001);
      order_parameters = ising_2d([temperature], 'chi', chi_values(c), 'N', N_values(n));
      dataset(c, n) = order_parameters(1);
    end
    plot(chi_values, dataset(:, n), MARKERS(n));
  end

  xlabel('\chi');
  ylabel('|m|');
  title('Order parameter at T_c')
  legend(legend_labels);
  hold off
end
