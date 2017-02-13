function [C, T] = calculate_environment_fixed_N(temperature, chi, N, initial_C, initial_T)
  C = initial_C;
  T = initial_T;

  for iteration = 1:N
    % [C, T, singular_values] = grow_lattice(C, T, a, chi);
    [C, T, singular_values] = grow_lattice(temperature, chi, C, T);
  end
end
