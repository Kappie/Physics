function set_markers(handles)
  MARKERS = markers();
  for i = 1:numel(handles)
    set(handles(i), 'marker', MARKERS(i), 'LineStyle', 'none')
  end
end
