function m = exact_order_parameter(beta)
  beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
  if beta > beta_crit
    m = (1 - sinh(2*beta).^-4).^(1/8);
  else
    m = 0;
  end
end
