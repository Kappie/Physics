function points = rounded_powerspace(power, startpoint, endpoint, number_of_points)
  % Analogous to linspace, creates a vector that is equally spaced when raised
  % to a power.

  points = fliplr( arrayfun(@round, linspace(startpoint, endpoint, number_of_points) .^ (1/power)) );

end
