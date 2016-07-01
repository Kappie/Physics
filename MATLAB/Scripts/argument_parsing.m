function argument_parsing(temperatures, varargin)
  p = inputParser;
  % default_chi = 4;
  % default_tolerance = 1e-6;

  addRequired(p, 'temperatures');
  addParameter(p, 'chi', 4);
  addParameter(p, 'tolerance', 1e-6);

  parse(p, temperatures, varargin{:});

  chi = p.Results.chi
  tolerance = p.Results.tolerance
  temperatures


end
