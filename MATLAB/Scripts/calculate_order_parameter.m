function data_points = calculate_order_parameter( temperatures, chi_values, N_values )
  data_points = zeros(numel(temperatures), numel(chi_values), numel(N_values));
  environment_tensors = find_or_calculate_environment_tensors(temperatures, chi_values, N_values);

  for t = 1:numel(temperatures)
    for c = 1:numel(chi_values)
      for n = 1:numel(N_values)
        C = environment_tensors(t, c, n).C;
        T = environment_tensors(t, c, n).T;
        data_points(t, c, n) = order_parameter(temperatures(t), C, T);
      end
    end
  end
end
