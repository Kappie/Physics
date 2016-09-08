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
end

function tensor = find_or_calculate_environment(temperature, chi, N)
  % Returns a struct with fields C and T containing the converged environment tensors.
  save_to_db = true;
  tensor = struct('C', [], 'T', []);

  initial_C = spin_up_initial_C(temperature);
  initial_T = spin_up_initial_T(temperature);
  number_of_iterations_remaining = N;

  query = ['SELECT * ' ...
    'FROM tensors ' ...
    'WHERE temperature = ? AND chi = ? AND n <= ? ' ...
    'ORDER BY n DESC ' ...
    'LIMIT 1'];
  found_record = sqlite3.execute(query, temperature, chi, N);

  if ~isempty(found_record)
    [initial_C, initial_T] = deserialize_tensors(found_record);
    number_of_iterations_remaining = number_of_iterations_remaining - found_record.n;
    display('Found a record:')
    display(['number of iterations remaining: ' num2str(number_of_iterations_remaining)]);
  end

  if number_of_iterations_remaining == 0
    tensor.C = initial_C; tensor.T = initial_T;
  else
    [tensor.C, tensor.T] = calculate_environment(temperature, chi, N, initial_C, initial_T);
    if save_to_db
      sqlite3.execute('INSERT INTO tensors (temperature, chi, n, c, t) VALUES (?, ?, ?, ?, ?)', ...
        temperature, chi, N, getByteStreamFromArray(tensor.C), getByteStreamFromArray(tensor.T));
      display('saved stuff to db')
    end
  end
end
