function make_legend(values, variable_name)
  legend_labels = arrayfun(@(value) [variable_name ' = ' num2str(value)], values, 'UniformOutput', false);
  legend(legend_labels, 'Location', 'northwest');
end
