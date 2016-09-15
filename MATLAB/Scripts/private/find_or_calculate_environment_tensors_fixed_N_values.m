function tensors = find_or_calculate_environment_tensors_fixed_N_values(temperatures, chi_values, N_values)
  % Returns a struct array with structs with fields C and T containing the converged environment.
  sqlite3.open(DATABASE);

  for t = 1:numel(temperatures)
    for c = 1:numel(chi_values)
      for n  = 1:numel(N_values)
        tensors(t, c, n) = find_or_calculate_environment_fixed_N(temperatures(t), chi_values(c), N_values(n));
      end
    end
  end

  sqlite3.close(DATABASE);
end
