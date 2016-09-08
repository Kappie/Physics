function data_collapse
  % temperature_width = 0.005;
  % temperatures = linspace(T_crit - temperature_width, T_crit + temperature_width, 9);
  % N_values = 25:25:1000;
  % chi = 32;
  %
  % data_points = zeros(numel(temperatures), numel(N_values));
  %
  % for n = 1:numel(N_values)
  %   order_parameters = ising_2d(temperatures, 'chi', chi, 'N', N_values(n));
  %   data_points(:, n) = order_parameters;
  % end

  load('data_collapse.mat', 'temperatures', 'N_values', 'chi', 'data_points');
  MARKERS = markers();

  % Critical exponents
  % nu = 1, but we leave it out altogether.
  beta = 1/8;

  figure
  hold on

  for n = 1:numel(N_values)
    x_values = zeros(1, numel(temperatures));
    scaling_function_values = zeros(1, numel(temperatures));

    for t = 1:numel(temperatures)
      x_values(t) = reduced_T(temperatures(t)) * N_values(n);
      scaling_function_values(t) = N_values(n)^beta * data_points(t, n);
    end

    plot(x_values, scaling_function_values, MARKERS(mod(n, numel(MARKERS)) + 1))
  end

  xlabel('$ tN^{1/\nu} $', 'interpreter', 'latex')
  ylabel('$N^{\beta/\nu}m(t, N)$', 'interpreter', 'latex')
  title('Data collapse for scaling ansatz in N (number of steps). \chi = 32. 25 \leq N \leq 1000. |T - T_c| \leq 0.005')

  % make_legend(N_values, 'N')

  export_fig('../Plots/datacollapse_chi32_N25-1000.pdf')

  function t = reduced_T(T)
    t = (T - T_crit) / T_crit;
  end
end
