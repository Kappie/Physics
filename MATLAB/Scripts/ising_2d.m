function ising_2d
  % rng(10)

  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  J = 1;
  chi = 5;
  chi_init = 3;
  max_iterations = 200;
  tolerance = 1e-6;

  beta = 0.45;

  [C, T, iterations] = calculate_environment(beta, tolerance, max_iterations);

  m = order_parameter(beta, C, T)
  m_exact = exact_order_parameter(beta)

  % compare_environments(beta, C, T);

  % Works only for beta > beta_crit.
  function m = exact_order_parameter(beta)
    m = (1 - sinh(2*beta)^-4)^(1/8);
  end

  function m = order_parameter(beta, C, T)
    m = abs(magnetization(beta, C, T));
  end

  % Compute magnetization using converged environment tensors.
  function m = magnetization(beta, C, T)
    Z = partition_function(beta, C, T);
    unnormalized_magnetization = attach_environment(construct_b(beta), C, T);
    m = unnormalized_magnetization / Z;
  end

  function Z = partition_function(beta, C, T)
    Z = attach_environment(construct_a(beta), C, T);
  end

  function result = compare_environments(beta, C, T)
    a = construct_a(beta);
    z1 = attach_environment(a, C, T)

    env = environment(C, T);
    z2 = ncon({a, env}, {[1 2 3 4], [1 2 3 4]})

    result = 1;
  end

  % Expects either a or b.
  function result = attach_environment(tensor, C, T)
    env = environment(C, T);
    result = ncon({tensor, env}, {[1 2 3 4], [1 2 3 4]});
    %%%
    %%% Different way of calculating environment: gives the same result as current method.
    %%%
    % result = ncon({tensor, T, T, T, T}, ...
    %   {[1, 2, 3, 4], [1, -1, -2], [2, -3, -4], [3, -5, -6], [4, -7, -8]});
    % % Attach C.
    % result = ncon({result, C, C, C, C}, ...
    %   {[1, 2, 3, 4, 5, 6, 7, 8], [2, 3], [4, 5], [6, 7], [8, 1]});
  end

  function environment = environment(C, T)
    % Final order is such that the physical dimension of a quarter of the total
    % environment comes first. Second comes the T leg, third the C leg.
    quarter = ncon({C, T}, {[1, -1], [-2, 1, -3]}, [1], [-2 -3 -1]);
    % Final order is such that two physical dimensions come first, second the T leg,
    % third the C leg.
    half = ncon({quarter, quarter}, {[-1 1 -2], [-3 -4 1]}, [1], [-1 -3 -4 -2]);
    environment = ncon({half, half}, {[-1 -2 1 2], [-3 -4 2 1]});
  end

  function [C, T, iteration] = calculate_environment(beta, tolerance, max_iterations)
    C = random_C();
    T = random_T();
    singular_values = initial_singular_values();
    a = construct_a(beta);

    for iteration = 1:max_iterations
      singular_values_old = singular_values;
      [C, T, singular_values] = grow_lattice(C, T, a);

      c = convergence(singular_values, singular_values_old);
      disp(iteration)
      disp(c)

      if convergence(singular_values, singular_values_old) < tolerance
        break
      end
    end
  end

  function [C, T, singular_values] = grow_lattice(C, T, a)
    % Final order is specified so that the new tensor is ordered according to
    % [d, chi, c, chi], with the pairs of c, chi corresponding to what will be the reshaped
    % legs of the new C.
    C = ncon({C, T, T, a}, {[1, 2], [3, 1, -1], [4, 2, -2], [3, -3, -4, 4]}, ...
    [1, 2, 3, 4], [-3 -1 -4 -2]);
    [U, s, U_transpose] = tensorsvd(C, [1 2], [3 4], chi, 'n');
    % Here, we only take the chi most relevant eigenvectors.
    C = ncon({C, U_transpose, U}, {[1 2 3 4], [1 2 -1], [3 4 -2]});

    % Again, final order is chosen such that we have [d, chi, d, chi, d], i.e.
    % left legs, right legs, middle physical leg.
    T = ncon({T, a}, {[1, -1, -2], [1, -3, -4, -5]}, [1], [-3 -1 -4 -2 -5]);
    % Again, keeping only chi most relevant eigenvectors, and being careful
    % to attach U first, so that U.U_transpose becomes a unity in the relevant
    % subspace when we construct the lattice later.
    T = ncon({T, U, U_transpose}, {[1 2 3 4 -1], [1 2 -2], [3 4 -3]});

    % Scale elements to prevent values from diverging when performing numerous growth steps.
    % Resymmetrize to prevent numerical errors adding up to unsymmetrize tensors.
    C = scale_by_largest_element(C);
    C = symmetrize_C(C);
    T = scale_by_largest_element(T);
    T = symmetrize_T(T);
    singular_values = scale_singular_values(diag(s));
  end

  function c = convergence(singular_values, singular_values_old)
    % Sometimes it happens that the current singular values vector is smaller
    % than the old one, because MATLAB's svd procedure throws away excessively
    % small singular values. The line below adds zeros to singular_values to match
    % the dimension of singular_values_old.
    singular_values(numel(singular_values_old)) = 0;
    c = sum(abs(singular_values - singular_values_old));
  end

  function M = scale_by_largest_element(M)
    M = M / max(M(:));
  end

  function s = scale_singular_values(singular_values)
    s = singular_values / sum(singular_values);
  end

  function Q = construct_Q(beta)
    Q = [exp(beta*J) exp(-beta*J); exp(-beta*J) exp(beta*J)];
  end

  function a = construct_a(beta)
    delta = construct_delta();
    % We need square root of a matrix here, not the square root of the elements!
    P = sqrtm(construct_Q(beta));
    a = ncon({P, P, P, P, delta}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
  end

  function b = construct_b(beta)
    g = construct_g();
    P = sqrtm(construct_Q(beta));
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

  function s = initial_singular_values()
    s = ones(chi, 1) / chi;
  end

  function T = random_T()
    T = symmetrize_T(rand(2, chi_init, chi_init));
  end

  function C = random_C()
    C = symmetrize_C(rand(chi_init));
  end

  function C = symmetrize_C(C)
    C = symmetrize(C);
  end

  function T = symmetrize_T(T)
    % Why do we not symmetrize in the physical dimension?
    % Squeeze deletes the singleton dimension to obtain a chi x chi matrix.
    T(1,:,:) = symmetrize(squeeze(T(1,:,:)));
    T(2,:,:) = symmetrize(squeeze(T(2,:,:)));
  end

  function m = symmetrize(m)
    m = triu(m) + triu(m, 1)';
  end
end
