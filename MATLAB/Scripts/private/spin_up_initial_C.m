function C = spin_up_initial_C(temperature)
  spin_up_tensor = corner_delta();
  spin_up_tensor(2, 2) = 0;
  P = construct_P(temperature);
  C = ncon({P, P, spin_up_tensor}, {[-1 1], [-2 2], [1 2]});
end
