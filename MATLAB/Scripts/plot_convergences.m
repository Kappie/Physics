beta_crit = log(1 + sqrt(2)) / 2; % ~0.44
T_crit = 1 / beta_crit;

data_dir = '~/Documents/Natuurkunde/Scriptie/Code/Data/2D_Ising/convergences/';
file_names = { ...
  'convergences_chi16T2.1692.dat', ...
  'convergences_chi16T2.2192.dat', ...
  'convergences_chi16T2.2692.dat', ...
  'convergences_chi16T2.3192.dat', ...
  'convergences_chi16T2.3692.dat' };



figure

xlabel('iterations');
ylabel('convergence');

for name = file_names
  full_name = fullfile(data_dir, name);
  data = dlmread(full_name{1});
  semilogy(data)
  hold on
end

hold off;
