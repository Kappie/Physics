function P = construct_P(temperature)
  % We need square root of a matrix here, not the square root of the elements!
  P = sqrtm(construct_Q(temperature));
end

function Q = construct_Q(temperature)
  Q = [exp((1/temperature)*J) exp(-(1/temperature)*J); exp(-(1/temperature)*J) exp((1/temperature)*J)];
end
