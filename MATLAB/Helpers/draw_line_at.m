function draw_line_at(T)
  % Yeah this is gonna introduce some bug at some point.
  axis manual
  line([T, T], [0, 1], 'LineStyle', '--');
end
