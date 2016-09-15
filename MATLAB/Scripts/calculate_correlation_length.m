function data_points = calculate_correlation_length( temperatures, chi_values, varargin )
  data_points = calculate_quantity( @correlation_length, temperatures, chi_values, varargin );
end
