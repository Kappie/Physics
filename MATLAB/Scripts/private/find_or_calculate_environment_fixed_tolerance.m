function tensors = find_or_calculate_environment_fixed_tolerance(temperature, chi, tolerance)
  % Returns a struct with fields C and T containing the converged environment tensors.
  tensors = struct('C', [], 'T', []);

  initial_C = spin_up_initial_C(temperature);
  initial_T = spin_up_initial_T(temperature);

  % I look for all records with the same temperature, lesser or equal chi and greater or equal tolerance.
  % If I find an exact match (same temperature, chi, tolerance as I'm trying to simulate)
  % I do not simulate again and just return the C, T tensors from the database.
  % If I find a record with matching temperature and lesser chi or higher tolerance (highest chi takes precedence)
  % I select the C, T from that record to use as initial C, T for the new simulation.
  simulation = true;
  query = ['SELECT * ' ...
    'FROM tensors ' ...
    'WHERE temperature = ? AND chi <= ? AND tolerance >= ? ' ...
    'ORDER BY chi DESC, tolerance ASC ' ...
    'LIMIT 1'];
  found_record = sqlite3.execute(query, temperature, chi, tolerance);

  if ~isempty(found_record)
    [found_C, found_T] = deserialize_tensors(found_record);

    % Found exact tensors I was looking for; do not simulate at all.
    if found_record.chi == chi & found_record.tolerance == tolerance
      simulation = false;
      tensors.C = found_C;
      tensors.T = found_T;
      display('I loaded stuff from the DB.')
    % Found record with higher tolerance. Simulate with found tensors as initial conditions.
    else
      initial_C = found_C;
      initial_T = found_T;
    end
  end

  if simulation
    [tensors.C, tensors.T, converged] = calculate_environment_fixed_tolerance(temperature, chi, tolerance, initial_C, initial_T);

    if converged && SAVE_TO_DB
      serialized_C = getByteStreamFromArray(tensors.C);
      serialized_T = getByteStreamFromArray(tensors.T);
      sqlite3.execute('INSERT INTO tensors (c, t, temperature, chi, tolerance) VALUES (?, ?, ?, ?, ?)', ...
        serialized_C, serialized_T, temperature, chi, tolerance);
      display('I put stuff in the DB:')
    elseif ~converged
      display('Failed to converge:')
    elseif ~SAVE_TO_DB
      display('Not saving to db.')
    end

    display(['temp = ' num2str(temperature) ' chi = ' num2str(chi) ' tolerance = ' num2str(tolerance)])
  end
end
