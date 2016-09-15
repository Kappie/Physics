function [C, T, converged] = calculate_environment_fixed_tolerance(temperature, chi, tolerance, initial_C, initial_T)
  MAX_ITERATIONS = 1e5;
  C = initial_C;
  T = initial_T;
  singular_values = initial_singular_values(chi);
  converged = false;

  for iteration = 1:MAX_ITERATIONS
    singular_values_old = singular_values;
    [C, T, singular_values] = grow_lattice(temperature, chi, C, T);

    if convergence(singular_values, singular_values_old, chi) < tolerance
      converged = true;
      break;
    end
  end
end

function s = initial_singular_values(chi)
  s = ones(chi, 1) / chi;
end

function c = convergence(singular_values, singular_values_old, chi)
  % TODO: something weird happened last time: when calculating correlation length
  % for T >> Tc, I found that the T tensor did not have the right dimension.

  % Sometimes it happens that the current singular values vector is smaller
  % than the old one, because MATLAB's svd procedure throws away excessively
  % small singular values. The code below adds zeros to singular_values to match
  % the dimension of singular_values_old.
  if size(singular_values, 1) < chi
    singular_values(chi) = 0;
  end

  % If chi_init is small enough, the bond dimension of C and T will not exceed
  % chi for the first few steps.
  if size(singular_values_old, 1) < chi
    singular_values_old(chi) = 0;
  end

  c = sum(abs(singular_values - singular_values_old));
end
