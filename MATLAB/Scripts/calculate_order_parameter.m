function data_points = calculate_order_parameter( temperatures, chi_values, N_values )
  data_points = calculate_quantity( @order_parameter, temperatures, chi_values, N_values )
end
