function tensors = find_or_calculate_environment_fixed_N(temperature, chi, N)
  % Returns a struct with fields C and T containing the converged environment tensors.
  tensors = struct('C', [], 'T', []);

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
  else
    display('Did not find a record')
  end

  if number_of_iterations_remaining == 0
    tensors.C = initial_C; tensors.T = initial_T;
  else
    [tensors.C, tensors.T] = calculate_environment_fixed_N(temperature, chi, number_of_iterations_remaining, initial_C, initial_T);
    if SAVE_TO_DB
      sqlite3.execute('INSERT INTO tensors (temperature, chi, n, c, t) VALUES (?, ?, ?, ?, ?)', ...
        temperature, chi, N, getByteStreamFromArray(tensors.C), getByteStreamFromArray(tensors.T));
      display('saved stuff to db')
    end
  end
end
