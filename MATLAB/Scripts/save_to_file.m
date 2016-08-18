function save_to_file(columns, path, overwrite)
  file_exists = exist(path, 'file');
  if overwrite || ~file_exists
    dlmwrite(path, columns, 'delimiter', ' ', 'precision', 10);
    sprintf(['Wrote to file ' path '.\n']);
  else
    sprintf('Not overwriting, file exists.');
  end
end
