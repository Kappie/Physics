function correlation_length = correlation_length(temperature, chi, N)
  environment = find_or_calculate_environment(temperature, chi, N);
  T = environment.T;
  transfer_matrix = ncon({T, T}, {[1 -1 -2], [1 -3 -4]});
  % reshape into 2 * chi x 2 * chi matrix
  [transfer_matrix, ~, ~] = lreshape(transfer_matrix, [1 2], [3 4]);
  eigenvalues = eigs(transfer_matrix, 2);
  correlation_length = 1 / log(eigenvalues(1) / eigenvalues(2));
end
