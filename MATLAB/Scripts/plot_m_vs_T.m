function plot_m_vs_T

  % tolerance = 1e-7;
  % max_iterations = 1e9;

  N = 1000;

  width = 0.005;
  temperatures = linspace(T_crit - width, T_crit + width, 9);
  chi_values = [32];
  % chi_values = [2, 4];

  legend_labels = arrayfun(@(chi) ['chi = ' num2str(chi)], chi_values, 'UniformOutput', false);
  MARKERS = markers();

  figure
  hold on

  for i = 1:numel(chi_values)
    % order_parameters = ising_2d(temperatures, 'chi', chi_values(i), ...
    %   'tolerance', tolerance, 'max_iterations', max_iterations);
    order_parameters = ising_2d(temperatures, 'chi', chi_values(i), ...
      'N', N);
    plot(temperatures, order_parameters, ['--', MARKERS(i)])
  end

  legend(legend_labels, 'Location', 'southwest')
  line([T_crit, T_crit], [0, 1], 'LineStyle', '--');
  xlabel('T')
  ylabel('|m|')
  title('Magnetization of 2D Ising model at tolerance = 1e-7')
  hold off

  export_fig('test.pdf')
end
