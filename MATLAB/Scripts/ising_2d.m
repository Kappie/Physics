function ising_2d

  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  J = 1;
  chi = 5;
  chi_init = 2;
  construct_a(0.35)

  % Works only for beta > beta_crit.
  function m = exact_magnetization(beta)
    m = (1 - sinh(2*beta)^-4)^(1/8);
  end

  function [T C] = perform_ctm_step()
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
    g = construct_g()
    P = sqrt(construct_Q(beta));
    b = ncon({P, P, P, P, g}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]})
  end

  function delta = construct_delta()
    delta = zeros(2, 2, 2, 2);
    delta(1, 1, 1, 1) = 1;
    delta(2, 2, 2, 2) = 1;
  end

  function g = construct_g()
    g = construct_delta()
    g(2, 2, 2, 2) = -1
  end

  function t = random_t
  end

  function c = random_c
  end
end
