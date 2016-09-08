function tensors = find_or_calculate_environment_tensors(temperatures, chi_values, N_values)
  % Returns a struct array with fields C and T containing the converged environment.
  % I look for records with the same temperature, same chi, and lower N.
  % If I find a such a tensor, I subtract the number of steps I need to take in my simulation.

  % tensors(1:numel(temperatures), 1:numel(chi_values), 1:numel(N_values)) = struct('C', {}, 'T', {});

  DATABASE = 'converged_tensors.db';
  sqlite3.open(DATABASE);

  for t = 1:numel(temperatures)
    for c = 1:numel(chi_values)
      for n  = 1:numel(N_values)
        tensors(t, c, n) = find_or_calculate_environment(temperatures(t), chi_values(c), N_values(n));
      end
    end
  end

  sqlite3.close(DATABASE);
end
