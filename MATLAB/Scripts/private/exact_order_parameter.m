function m = exact_order_parameter(temperature)
  if temperature < T_crit
    m = (1 - sinh(2*beta).^-4).^(1/8);
  else
    m = 0;
  end
end
