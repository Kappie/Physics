function data_points = calculate_order_parameter( temperatures, chi_values, varargin )
  data_points = calculate_quantity( @order_parameter, temperatures, chi_values, varargin );
end
