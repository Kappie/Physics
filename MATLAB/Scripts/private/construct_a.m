function a = construct_a(temperature)
  delta = construct_delta();
  P = construct_P(temperature);
  a = ncon({P, P, P, P, delta}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
end
