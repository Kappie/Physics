function ising_2d
beta_crit = log(1 + sqrt(2)) / 2;

end

function m = exact_magnetization(beta)
m = (1 - sinh(2*beta)^-4)^(1/8);
end

function Q = Q(s1, s2, beta, J)
Q = exp(beta * J * s1 * s2);
end