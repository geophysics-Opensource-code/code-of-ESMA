function y = calculate_diffusive_wave_Hermite_adaptive(x, t, C, cgamma, z_asym, vareps)
%% calculate_diffusive_wave_solution: calculate the exact solution of the diffusive wave equation

% Description:
   % [xmin, xmax]: the range of spatial domain
   % Nx:           the size of collocation points
   % tfin:         the final time
   % C:            the group velocity
   % cgamma:       the order of the Caputo fractional derivate
   % MxG:          the points of Laguerre-Gauss nodes
   % z_asym :      the threshold to truncate the Green function

%% Grid mesh

tic
%% Preparing the Mainardi function
s = sqrt(C)*t^(1-cgamma);
y1 = -1000;
y2 = 0;
MxH = 100;
vareps2 = 10^-5;

while (abs(y1 - y2) > vareps2) 
  y1 = y2;
  MxH = MxH + 25;
  hermite_rule ( MxH, 0, 1, 0, 'HGQ' );
  Y  = load('HGQ_x.txt');
  WT = load('HGQ_w.txt');
  delete('HGQ_*')

  %% Calculation
  y2 = 0;
  for j = 1:MxH
    MF = calculate_Mainardi_function((x - Y(j))/s, 1-cgamma, z_asym, vareps);
    y2 = y2 + WT(j)*MF;
  end
  y2 = y2/(2*s);
  disp([num2str(MxH), ' Value: ', num2str(y2, '% 20.16f')])
end
toc

y = y2;

end