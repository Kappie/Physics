function f = exact_free_energy_per_site(temperature)
  J = 1;
  k = sinh(2*(1/temperature)*J)^-2;

  function int = integrand(theta)
    int = log( cosh(2*(1/temperature)*J)^2 + (1/k)*sqrt(1+k^2-2*k*cos(2*theta)) );
  end

  f = (-1/(1/temperature))*( log(2)/2 + (1/(2*pi))*integral(@integrand, 0, pi) );
end
