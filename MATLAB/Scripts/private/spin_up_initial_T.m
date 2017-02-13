function T = spin_up_initial_T(temperature)
  spin_up_tensor = edge_delta();
  spin_up_tensor(2, 2, 2) = 0;
  P = construct_P(temperature);
  T = ncon({P, P, P, spin_up_tensor}, {[-1 1], [-2 2], [-3 3], [1 2 3]});
end
