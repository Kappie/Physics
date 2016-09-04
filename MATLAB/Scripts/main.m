function main
  clc;
  format long;

  % plot_m_vs_chi_at_tcrit;
% plot_m_vs_N_at_tcrit
  plot_free_energy_vs_N_at_tcrit
  % plot_m_vs_T


  % max_iterations = 1e6;
  % tolerance = 1e-7;

  % result = ising_2d([temperature], 'chi', chi, 'N', N)

  % EXPERIMENT

  % EXPERIMENT
  % Plot magnetization versus 1/chi for high chi_values.

  % chi_values = [2, 4, 8, 16, 32, 64, 80, 96, 112];
  % tolerance_values = [1e-4, 1e-5, 1e-6, 1e-7];
  % temperature = T_crit
  % legend_labels = arrayfun(@(tolerance) ['tolerance = ' num2str(tolerance, '%.0e')], tolerance_values, 'UniformOutput', false);
  %
  % magnetizations = zeros(numel(chi_values), numel(tolerance_values));
  %
  % for t = 1:numel(tolerance_values)
  %   for c = 1:numel(chi_values)
  %     profile on
  %     result = ising_2d([temperature], 'chi', chi_values(c), 'tolerance', ...
  %       tolerance_values(t), 'max_iterations', max_iterations);
  %     profsave
  %     magnetizations(c, t) = result(1);
  %   end
  % end

  % figure
  % plot(1./chi_values, magnetizations, 'o--')
  % legend(legend_labels);
  % xlabel('1 / chi')
  % ylabel('|m|')
  % title('Magnetization of 2D Ising model at T_{crit}')

end
