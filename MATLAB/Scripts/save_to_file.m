function save_to_file(columns, path)
  write = true;
  if exist(path, 'file')
    answer = input('Warning, file exists. Do you want to overwrite? y/n [n]\n', 's');
    if ~strcmp(answer, 'y')
      write = false;
    end
  end

  if write
    dlmwrite(path, columns, 'delimiter', ' ', 'precision', 5);
    sprintf(['Wrote to file ' path '.\n']);
  end
end
