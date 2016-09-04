function plot_free_energy_vs_N_at_tcrit
  temperature = T_crit;
  chi_values = [4, 16, 32];
  N_values = 50:50:1000;
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
    semilogy(N_values, dataset(:, c), MARKERS(c))
    hold on
  end

  legend(legend_labels)
  xlabel('N')
  ylabel('f')
  title('Free energy per site at T_c')
  hold off
end
