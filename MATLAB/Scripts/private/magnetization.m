function m = magnetization(temperature, C, T)
  Z = partition_function(temperature, C, T);
  unnormalized_magnetization = attach_environment(construct_b(temperature), C, T);
  m = unnormalized_magnetization / Z;
end
