function plot_free_energy_vs_chi
  temperature = T_crit;
  chi_values = [4:32];
  N_values = [500, 7500, 12500];
  exact_value = exact_free_energy_per_site(1/temperature);

  legend_labels = arrayfun(@(N) ['N = ' num2str(N)], N_values, 'UniformOutput', false);
  MARKERS = markers;

  dataset = zeros(numel(chi_values), numel(N_values));

  for n = 1:numel(N_values)
    for c = 1:numel(chi_values)
      free_energies = ising_2d([temperature], 'chi', chi_values(c), 'N', N_values(n));
      dataset(c, n) = abs(free_energies(1) - exact_value);
    end
    loglog(chi_values, dataset(:, n), MARKERS(n));
    % save('free_energies_chi_values_N7500.mat', 'chi_values', 'dataset')
    hold on
  end

  xlabel('\chi');
  ylabel('$| f(N, \chi) - f(\infty) |$', 'interpreter', 'latex');
  title('Error in free energy per site at T_{crit}')
  legend(legend_labels);
  hold off
end
