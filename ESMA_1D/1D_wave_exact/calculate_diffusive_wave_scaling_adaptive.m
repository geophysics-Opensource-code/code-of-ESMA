function y = calculate_diffusive_wave_scaling_adaptive(x, t, C, cgamma, x_asym, vareps)
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
MxG = 100;
vareps2 = 10^-6;

while (abs(y1 - y2) > vareps2) 
  y1 = y2;
  MxG = MxG + 25;
  [Y, WT] = Laguerre_Gauss_FSU( MxG, 0, 4); 
  MF = zeros(1, MxG);
  for j = 1:MxG
    MF(j) = calculate_Mainardi_function(Y(j)/s, 1-cgamma, x_asym, vareps);
  end
  
  %% calculation
  y2 = 0;
  for j = 1:MxG
    z = (exp(-(x-Y(j))^2) + exp(-(x+Y(j))^2));
    y2 = y2 + z*WT(j)*MF(j);
  end
  y2 = y2/(2*s);
  disp([num2str(MxG), ' Value: ', num2str(y2, '% 20.16f')])
end
toc

y = y2;

end