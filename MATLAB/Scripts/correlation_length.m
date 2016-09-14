function correlation_length = correlation_length(temperature, ~, T)
  % construct row of spins, and attach two rows together.
  transfer_matrix = ncon({T, T}, {[1 -1 -3], [1 -2 -4]});
  transfer_matrix = ncon({transfer_matrix, transfer_matrix}, {[1 2 -1 -2], [1 2 -3 -4]});

  % reshape into chi^2 x chi^2 matrix
  [transfer_matrix, ~, ~] = lreshape(transfer_matrix, [1 2], [3 4]);

  eigenvalues = eigs(transfer_matrix, 2);
  correlation_length = 1 / log(eigenvalues(1) / eigenvalues(2));
end
