function strain = calculate_diffusive_wave_Hermite(xmesh, t, C, cgamma, MxH, z_asym, vareps)
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
Nx = size(xmesh, 2);
strain = zeros(1, Nx);

tic
%% Preparing the Mainardi function
s = sqrt(C)*t^(1-cgamma);

hermite_rule ( MxH, 0, 1, 0, 'HGQ' );
Y  = load('HGQ_x.txt');
WT = load('HGQ_w.txt');
delete('HGQ_*')

MF = zeros(1, MxH);
for j = 1:MxH
  MF(j) = calculate_Mainardi_function(Y(j), 1-cgamma, z_asym, vareps);
end

%% Calculation
for i = 1:Nx
  for j = 1:MxH
    MF = calculate_Mainardi_function((xmesh(i) - Y(j))/s, 1-cgamma, z_asym, vareps);
    strain(i) = strain(i) + WT(j)*MF;
  end
  disp([num2str(xmesh(i)), ' value: ', num2str(strain(i)/(2*s))])
end
strain = strain/(2*s);
disp(['Current time: ', num2str(t), ' Positive Peak: ', num2str(max(strain), '% 20.16f')])
toc

%% Plot
% figure 
% box on 
% plot(xmesh, strain, '-o', 'LineWidth', 1, 'Markersize', 6);
% xlabel('$$x$$', 'fontsize', 20, 'Interpreter', 'latex')
% ylabel('strain', 'fontsize', 20, 'Interpreter', 'latex')
% axis([xmesh(1) xmesh(end) min(min(strain)) max(max(strain))])
% drawnow;

% %% Save
% save varepsilon_cgamma_t varepsilon xmesh 
end