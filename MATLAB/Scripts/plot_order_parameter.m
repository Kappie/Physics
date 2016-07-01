function plot_order_parameter(temperatures, order_parameters)
  plt = Plot(temperatures, order_parameters);
  plt.LineStyle = {'none'};
  % plt.Markers = {'o'};

  plt.Title = 'Order parameter';
  plt.XLabel = 'T';
  plt.YLabel = 'm';
  plt.Markers = {'d'};
  plt.ShowBox = 'off';
  plt.LineStyle = '--';

  % f = 50;  % frequency
  % Vm = 10; % peak
  % phi = 0; % phase
  %
  % % generate the signal
  % t = [0:0.0001:3/f];
  % th = 2*pi*f*t;
  % v = Vm*sin(th+phi);
  %
  % % plot it
  % plt = Plot(t*1E3, v);
  %
  % plt.Title = 'Voltage as a function of time'; % plot title
  % plt.XLabel = 'Time, t (ms)'; % xlabel
  % plt.YLabel = 'Voltage, V (V)'; %ylabel

  % plt.LineStyle = {'--'};
  % plt.export('plotSimple1.png');
end
