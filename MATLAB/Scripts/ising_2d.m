function result = ising_2d(temperatures, varargin)
  p = inputParser;
  default_chi = [];
  default_tolerance = [];
  default_N = [];

  default_chi_init = 2;
  default_max_iterations = 1e6;
  default_min_iterations = 0;
  default_tensor_initialization = 'spin-up';
  default_traversal_order = 'standard';

  addRequired(p, 'temperatures');
  addParameter(p, 'chi', default_chi);
  addParameter(p, 'chi_init', default_chi_init);
  addParameter(p, 'tolerance', default_tolerance);
  addParameter(p, 'max_iterations', default_max_iterations);
  addParameter(p, 'min_iterations', default_min_iterations);
  addParameter(p, 'tensor_initialization', default_tensor_initialization);
  addParameter(p, 'traversal_order', default_traversal_order);
  addParameter(p, 'N', default_N);

  parse(p, temperatures, varargin{:});

  chi = p.Results.chi;
  tolerance = p.Results.tolerance;
  N = p.Results.N;

  chi_init = p.Results.chi_init;
  max_iterations = p.Results.max_iterations;
  min_iterations = p.Results.min_iterations;
  tensor_initialization = p.Results.tensor_initialization;
  traversal_order = p.Results.traversal_order;

  database = 'converged_tensors.db';
  save_to_db = true;

  J = 1;

  % result = calculate_order_parameters(1./temperatures);
  result = calculate_free_energy_per_site(1./temperatures);

  function order_parameters = calculate_order_parameters(betas)
    number_of_points = numel(betas);
    order_parameters = zeros(1, number_of_points);

    for i = 1:number_of_points
      [C, T] = calculate_environment_if_it_does_not_exist(betas(i));
      order_parameters(i) = order_parameter(betas(i), C, T);
    end
  end

  function data_points = calculate_free_energy_per_site(betas)
    number_of_points = numel(betas);
    data_points = zeros(1, number_of_points);

    for i = 1:number_of_points
      [C, T] = calculate_environment_if_it_does_not_exist(betas(i));
      data_points(i) = free_energy_per_site(betas(i), C, T);
    end
  end

  function [C, T] = calculate_environment_if_it_does_not_exist(beta)
    % This function tries to find existing converged environment tensors in a sqlite3 database.
    % If the tensors I need do not exist, I calculate them. I try various strategies (depending on
    % whether I do a fixed N simulation or try to reach a certain convergence, to come up with
    % initial tensors I can use from the database.
    % Database schema: CREATE TABLE tensors (c BLOB, t BLOB, temperature NUMERIC, chi NUMERIC, tolerance NUMERIC, n NUMERIC)
    initial_C = spin_up_initial_C(beta);
    initial_T = spin_up_initial_T(beta);


    % Without pause, database can crash if I make too many requests.
    % This is really just a poor hack because I didn't organize everything into one query yet.
    pause(0.001);
    sqlite3.open(database);

    if N
      [C, T] = environment_with_fixed_iterations(beta, initial_C, initial_T);
    else
      [C, T] = environment_with_fixed_tolerance(beta, initial_C, initial_T);
    end
  end

  function [C, T] = environment_with_fixed_iterations(beta, initial_C, initial_T)
    % I look for records with the same temperature, same chi, and lower N.
    % If I find a such a tensor, I subtract the number of steps I need to take in my simulation.
    query = ['SELECT * ' ...
      'FROM tensors ' ...
      'WHERE temperature = ? AND chi = ? AND n <= ? ' ...
      'ORDER BY n DESC ' ...
      'LIMIT 1'];
    found_record = sqlite3.execute(query, 1/beta, chi, N);

    number_of_iterations_remaining = N;

    if ~isempty(found_record)
      [initial_C, initial_T] = deserialize_tensors(found_record);
      number_of_iterations_remaining = number_of_iterations_remaining - found_record.n;
      display('Found a record:')
      display(['number of iterations remaining: ' num2str(number_of_iterations_remaining)]);
    end

    if number_of_iterations_remaining == 0
      C = initial_C; T = initial_T;
    else
      [C, T, ~] = calculate_environment(beta, initial_C, initial_T, number_of_iterations_remaining);
      if save_to_db
        sqlite3.execute('INSERT INTO tensors (temperature, chi, n, c, t) VALUES (?, ?, ?, ?, ?)', ...
          1/beta, chi, N, getByteStreamFromArray(C), getByteStreamFromArray(T));
        display('saved stuff to db')
      end
    end
  end

  function [C, T] = environment_with_fixed_tolerance(beta, initial_C, initial_T)
    % I look for all records with the same temperature, lesser or equal chi and greater or equal tolerance.
    % If I find an exact match (same temperature, chi, tolerance as I'm trying to simulate)
    % I do not simulate again and just return the C, T tensors from the database.
    % If I find a record with matching temperature and lesser chi or higher tolerance (highest chi takes precedence)
    % I select the C, T from that record to use as initial C, T for the new simulation.
    simulation = true;
    query = ['SELECT * ' ...
      'FROM tensors ' ...
      'WHERE temperature = ? AND chi <= ? AND tolerance >= ? ' ...
      'ORDER BY chi DESC, tolerance ASC ' ...
      'LIMIT 1'];
    found_record = sqlite3.execute(query, 1/beta, chi, tolerance);

    if ~isempty(found_record)
      [found_C, found_T] = deserialize_tensors(found_record);

      % Found exact tensors I was looking for; do not simulate at all.
      if found_record.chi == chi & found_record.tolerance == tolerance
        simulation = false;
        C = found_C;
        T = found_T;
        display('I loaded stuff from the DB.')
      else
        initial_C = found_C;
        initial_T = found_T;
      end
    end

    if simulation
      [C, T, converged] = calculate_environment(beta, initial_C, initial_T);

      if converged && save_to_db
        serialized_C = getByteStreamFromArray(C);
        serialized_T = getByteStreamFromArray(T);
        sqlite3.execute('INSERT INTO tensors (c, t, temperature, chi, tolerance) VALUES (?, ?, ?, ?, ?)', ...
          serialized_C, serialized_T, 1/beta, chi, tolerance);
        display('I put stuff in the DB:')
      elseif ~converged
        display('Failed to converge:')
      elseif ~save_to_db
        display('Not saving to db.')
      end

      display(['temp = ' num2str(1/beta) ' chi = ' num2str(chi) ' tolerance = ' num2str(tolerance)])
    end
  end

  function [C, T, converged] = calculate_environment(beta, initial_C, initial_T, varargin)
    p = inputParser;
    default_number_of_iterations = [];
    addRequired(p, 'beta');
    addRequired(p, 'initial_C');
    addRequired(p, 'initial_T');
    addOptional(p, 'number_of_iterations', default_number_of_iterations);

    parse(p, beta, initial_C, initial_T, varargin{:});
    number_of_iterations = p.Results.number_of_iterations;

    % Here, I use a hack to fix the number of iterations.
    if ~isequal(number_of_iterations, default_number_of_iterations)
      tolerance = -1;
      max_iterations = number_of_iterations;
      converged = [];
    end



    C = initial_C;
    T = initial_T;
    singular_values = initial_singular_values();
    a = construct_a(beta);

    for iteration = 1:max_iterations
      singular_values_old = singular_values;
      [C, T, singular_values] = grow_lattice(C, T, a, chi);
      c = convergence(singular_values, singular_values_old);

      if c < tolerance && iteration >= min_iterations
        converged = true;
        break
      end
    end

    if c > tolerance
      converged = false;
    end

    % figure;
    % xlabel('iterations')
    % title(['Convergence and order parameter at T = ' num2str(1/beta) ' chi = ' num2str(chi)])
    %
    % yyaxi s left;
    % semilogy(cell2mat(convergences));
    % ylabel('convergence')
    %
    % yyaxis right;
    % plot(cell2mat(order_parameters));
    % ylabel('|m|')

    % data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/convergences/';
    % file_name = ['convergences_chi' num2str(chi) 'T' num2str(1/beta) '.dat'];
    % name = fullfile(data_dir, file_name);
    % save_to_file(cell2mat(convergences)', name, true);
  end

  function m = order_parameter(beta, C, T)
    m = abs(magnetization(beta, C, T));
  end

  function m = magnetization(beta, C, T)
    Z = partition_function(beta, C, T);
    unnormalized_magnetization = attach_environment(construct_b(beta), C, T);
    m = unnormalized_magnetization / Z;
  end

  % How does this follow from definition <s_i s_j> ?
  function f = free_energy_per_site(beta, C, T)
    Z = partition_function(beta, C, T);
    four_corners = ncon({C, C, C, C}, {[1 2], [2 3], [3 4], [4 1]});
    two_rows = ncon({C, T, C, C, T, C}, {[1 2], [3 2 4], [4 5], [5 6], [3 6 7], [7 1]});
    kappa = Z * four_corners / two_rows^2;
    f = - (1/beta)*log(kappa);
  end

  function Z = partition_function(beta, C, T)
    Z = attach_environment(construct_a(beta), C, T);
  end

  % function result = compare_environments(beta, C, T)
  %   a = construct_a(beta);
  %   z1 = attach_environment(a, C, T)
  %
  %   env = environment(C, T);
  %   z2 = ncon({a, env}, {[1 2 3 4], [1 2 3 4]})
  %
  %   result = 1;
  % end

  % Expects either a or b.
  function result = attach_environment(tensor, C, T)
    env = environment(C, T);
    result = ncon({tensor, env}, {[1 2 3 4], [1 2 3 4]});
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



  function [C, T, singular_values] = grow_lattice(C, T, a, chi)
    % Final order is specified so that the new tensor is ordered according to
    % [d, chi, d, chi], with the pairs of d, chi corresponding to what will be the reshaped
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

  function M = scale_by_largest_element(M)
    M = M / max(M(:));
  end

  function s = scale_singular_values(singular_values)
    s = singular_values / sum(singular_values);
  end

  function Q = construct_Q(beta)
    Q = [exp(beta*J) exp(-beta*J); exp(-beta*J) exp(beta*J)];
  end

  function P = construct_P(beta)
    % We need square root of a matrix here, not the square root of the elements!
    P = sqrtm(construct_Q(beta));
  end

  function a = construct_a(beta)
    delta = construct_delta();
    P = construct_P(beta);
    a = ncon({P, P, P, P, delta}, {[-1, 1], [-2, 2], [-3, 3], [-4, 4], [1, 2, 3, 4]});
  end

  function b = construct_b(beta)
    g = construct_g();
    P = construct_P(beta);
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

  function delta = corner_delta()
    delta = zeros(2, 2);
    delta(1, 1) = 1;
    delta(2, 2) = 1;
  end

  function delta = edge_delta()
    delta = zeros(2, 2, 2);
    delta(1, 1, 1) = 1;
    delta(2, 2, 2) = 2;
  end

  function s = initial_singular_values()
    s = ones(chi, 1) / chi;
  end

  % Convention is: physical dimension first.
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

  function T = symmetric_initial_T(beta)
    delta = edge_delta();
    P = construct_P(beta);
    T = ncon({P, P, P, delta}, {[-1, 1], [-2, 2], [-3, 3], [1, 2, 3]});
  end

  function C = symmetric_initial_C(beta)
    delta = corner_delta();
    P = construct_P(beta);
    C = ncon({P, P, delta}, {[-1, 1], [-2, 2], [1, 2]});
  end

  % This corresponds to a corner with an upspin. Why exactly?
  function C = spin_up_initial_C(beta)
    spin_up_tensor = corner_delta();
    spin_up_tensor(2, 2) = 0;
    P = construct_P(beta);
    C = ncon({P, P, spin_up_tensor}, {[-1 1], [-2 2], [1 2]});
  end

  function T = spin_up_initial_T(beta)
    spin_up_tensor = edge_delta();
    spin_up_tensor(2, 2, 2) = 0;
    P = construct_P(beta);
    T = ncon({P, P, P, spin_up_tensor}, {[-1 1], [-2 2], [-3 3], [1 2 3]});
  end

  function [C, T] = deserialize_tensors(record)
    C = getArrayFromByteStream(record.c);
    T = getArrayFromByteStream(record.t);
  end
end
