function b = construct_b(temperature)
  g = construct_g();
  P = construct_P(temperature);
  b = ncon({P, P, P, P, g}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
end

function g = construct_g()
  g = construct_delta();
  g(2, 2, 2, 2) = -1;
end
