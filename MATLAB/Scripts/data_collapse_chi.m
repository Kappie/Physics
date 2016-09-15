function data_collapse_chi
  temperature_width = 0.001;
  T_pseudocrit = 2.270185314213022;
  temperatures = linspace(T_pseudocrit - temperature_width, T_pseudocrit + temperature_width, 11);
  tolerances = [1e-7];
  chi_values = [14 15 16];

  order_parameters = zeros(numel(temperatures), numel(chi_values));
  correlation_lengths = zeros(numel(temperatures), numel(chi_values));

  order_parameters = calculate_order_parameter(temperatures, chi_values, 'tolerances', tolerances);
  order_parameters = squeeze(order_parameters(:, :, 1));

  correlation_lengths = calculate_correlation_length(temperatures, chi_values, 'tolerances', tolerances);
  correlation_lengths = squeeze(correlation_lengths(:, :, 1));

  save('data_collapse_chi_zoom_t_pseudocrit.mat', 'temperatures', 'chi_values', 'order_parameters', 'correlation_lengths');

  load('data_collapse_chi_zoom_t_pseudocrit.mat', 'temperatures', 'chi_values', 'order_parameters', 'correlation_lengths');

  % Critical exponents
  % nu = 1, but we leave it out altogether.
  beta = 1/8;
  MARKERS = markers();

  figure
  hold on

  for c = 1:numel(chi_values)
    x_values = zeros(1, numel(temperatures));
    scaling_function_values = zeros(1, numel(temperatures));

    for t = 1:numel(temperatures)
      x_values(t) = reduced_T(temperatures(t)) * correlation_lengths(t, c);
      scaling_function_values(t) = order_parameters(t, c) * correlation_lengths(t, c)^beta;
    end

    plot(x_values, scaling_function_values, MARKERS(c));
  end

  make_legend(chi_values, 'chi')
  xlabel('$t\xi(\chi, t)^{1/\nu}$');
  ylabel('$m(t, \chi)\xi(\chi,t)^{\beta/\nu}$')

  function t = reduced_T(T)
    t = (T - T_crit) / T_crit;
  end
end
