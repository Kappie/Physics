function set_markers(handles)
  MARKERS = markers();
  for i = 1:numel(handles)
    marker_index = mod(i - 1, numel(MARKERS)) + 1;
    set(handles(i), 'marker', MARKERS(marker_index), 'LineStyle', 'none')
  end
end
