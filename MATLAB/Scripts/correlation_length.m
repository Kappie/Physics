function correlation_length = correlation_length(temperature, ~, T)
  % construct row of spins, and attach two rows together.
  transfer_matrix = ncon({T, construct_a(temperature), T}, {[1 -1 -4], [1 -2 2 -5], [2 -3 -6]});

  % reshape into 2chi^2 x 2chi^2 matrix
  [transfer_matrix, ~, ~] = lreshape(transfer_matrix, [1 2 3], [4 5 6]);

  eigenvalues = eigs(transfer_matrix, 2);
  correlation_length = 1 / log(eigenvalues(1) / eigenvalues(2));
end
