function data_points = calculate_correlation_length( temperatures, chi_values, N_values )
  data_points = calculate_quantity( @correlation_length, temperatures, chi_values, N_values );
end
