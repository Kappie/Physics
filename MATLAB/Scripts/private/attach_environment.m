function result = attach_environment(tensor, C, T)
  env = environment(C, T);
  result = ncon({tensor, env}, {[1 2 3 4], [1 2 3 4]});
end
