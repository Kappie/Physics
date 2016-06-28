function ising_2d

  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  J = 1;
  chi = 5;
  chi_init = 3;

  C = random_C();
  T = random_T();
  a = construct_a(0.35);

  grow_C(C, T, a);

  % Works only for beta > beta_crit.
  function m = exact_magnetization(beta)
    m = (1 - sinh(2*beta)^-4)^(1/8);
  end

  function [T C] = perform_ctm_step()
  end

  function C = grow_C(C, T, a)
    % Final order is specified so that the new tensor is ordered according to
    % [d, chi, c, chi], with the pairs of c, chi corresponding to the reshaped
    % legs of the new C.
    C = ncon({C, T, T, a}, {[1, 2], [3, 1, -1], [4, 2, -2], [3, -3, -4, 4]}, ...
    [1, 2, 3, 4], [-3 -1 -4 -2]);
    combined_dimension = size(C, 1) * size(C, 2);
    C = reshape(C, combined_dimension, combined_dimension)


  end

  % Keep only chi states of largest weight by performing singular value decomposition.
  % Assumes A is symmetric.
  function truncate_matrix(A, max_dimension)
    [U, s, V] = svd(A)
  end

  function T = grow_T(T, a)
  end

  function Q = construct_Q(beta)
    Q = [exp(beta*J) exp(-beta*J); exp(-beta*J) exp(beta*J)];
  end

  function a = construct_a(beta)
    delta = construct_delta();
    P = sqrt(construct_Q(beta));
    a = ncon({P, P, P, P, delta}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
  end

  function b = construct_b(beta)
    g = construct_g();
    P = sqrt(construct_Q(beta));
    b = ncon({P, P, P, P, g}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
  end

  function delta = construct_delta()
    delta = zeros(2, 2, 2, 2);
    delta(1, 1, 1, 1) = 1;
    delta(2, 2, 2, 2) = 1;
  end

  function g = construct_g()
    g = construct_delta();
    g(2, 2, 2, 2) = -1;
  end

  function T = random_T
    % Why do we not symmetrize in the physical dimension?
    T = rand(2, chi_init, chi_init);
    T(1,:,:) = random_symmetric_matrix(chi_init);
    T(2,:,:) = random_symmetric_matrix(chi_init);
  end

  function C = random_C
    C = random_symmetric_matrix(chi_init);
  end

  function r = random_symmetric_matrix(dim)
    r = rand(dim);
    r = triu(r) + triu(r, 1)';
  end

end
