function environment = environment(C, T)
  % Final order is such that the physical dimension of a quarter of the total
  % environment comes first. Second comes the T leg, third the C leg.
  quarter = ncon({C, T}, {[1, -1], [-2, 1, -3]}, [1], [-2 -3 -1]);
  % Final order is such that two physical dimensions come first, second the T leg,
  % third the C leg.
  half = ncon({quarter, quarter}, {[-1 1 -2], [-3 -4 1]}, [1], [-1 -3 -4 -2]);
  environment = ncon({half, half}, {[-1 -2 1 2], [-3 -4 2 1]});
end
