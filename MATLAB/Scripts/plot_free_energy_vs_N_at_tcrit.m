function plot_free_energy_vs_N_at_tcrit
  temperature = T_crit;
  chi_values = 4:4:32;
  N1 = 25:25:1000;
  N2 = arrayfun(@round, 10.^linspace(3, 4, 20));
  N_values = [N1'; N2'];
  exact_value = exact_free_energy_per_site(1/temperature);

  MARKERS = markers();
  legend_labels = arrayfun(@(chi) ['chi = ' num2str(chi)], chi_values, 'UniformOutput', false);

  dataset = zeros(numel(N_values), numel(chi_values));

  figure

  for c = 1:numel(chi_values)
    for n = 1:numel(N_values)
      free_energy_per_site = ising_2d([temperature], 'chi', chi_values(c), 'N', N_values(n));
      dataset(n,c) = abs(free_energy_per_site - exact_value);
    end
    % loglog(N_values, dataset(:, c), MARKERS(c))
    % hold on
  end

  save('free_energies_N_values_chi_values.mat', 'N_values', 'chi_values', 'dataset')

  % legend(legend_labels)
  % xlabel('N')
  % ylabel('$| f(N, \chi) - f(\infty) |$', 'interpreter', 'latex');
  % title('Error in free energy per site at T_{crit}')
  % hold off
end
