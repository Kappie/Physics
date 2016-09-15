function data_points = calculate_quantity( quantity, temperatures, chi_values, varargin )
  p = inputParser;
  addParameter(p, 'N_values', false);
  addParameter(p, 'tolerances', false);
  p.parse(varargin{:}{:});

  if p.Results.tolerances
    tolerances = p.Results.tolerances;
    data_points = calculate_quantity_fixed_tolerances(quantity, temperatures, chi_values, tolerances);
  elseif p.Results.N_values
    N_values = p.Results.N_values;
    data_points = calculate_quantity_fixed_N_values(quantity, temperatures, chi_values, N_values);
  else
    error('You did not specify either "tolerances" or "N_values" as a name-parameter pair.')
  end

end

function data_points = calculate_quantity_fixed_N_values(quantity, temperatures, chi_values, N_values)
  data_points = zeros(numel(temperatures), numel(chi_values), numel(N_values));
  environment_tensors = find_or_calculate_environment_tensors_fixed_N_values(temperatures, chi_values, N_values);

  for t = 1:numel(temperatures)
    for c = 1:numel(chi_values)
      for n = 1:numel(N_values)
        C = environment_tensors(t, c, n).C;
        T = environment_tensors(t, c, n).T;
        data_points(t, c, n) = quantity(temperatures(t), C, T);
      end
    end
  end
end

function data_points = calculate_quantity_fixed_tolerances(quantity, temperatures, chi_values, tolerances)
  data_points = zeros(numel(temperatures), numel(chi_values), numel(tolerances));
  environment_tensors = find_or_calculate_environment_tensors_fixed_tolerances(temperatures, chi_values, tolerances);

  for t = 1:numel(temperatures)
    for c = 1:numel(chi_values)
      for tol = 1:numel(tolerances)
        C = environment_tensors(t, c, tol).C;
        T = environment_tensors(t, c, tol).T;
        data_points(t, c, tol) = quantity(temperatures(t), C, T);
      end
    end
  end
end
