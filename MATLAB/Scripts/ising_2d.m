function ising_2d

J = -1;
beta_crit = log(1 + sqrt(2)) / 2;

end

function m = exact_magnetization(beta)
m = (1 - sinh(2*beta)^-4)^(1/8);
end