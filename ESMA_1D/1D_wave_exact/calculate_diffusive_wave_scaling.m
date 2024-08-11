function strain = calculate_diffusive_wave_scaling(xmesh, t, C, cgamma, MxG, z_asym, vareps)
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
[Y, WT] = Laguerre_Gauss_FSU( MxG, 0, 1/s);
MF = zeros(1, MxG);
for j = 1:MxG
  MF(j) = calculate_Mainardi_function(Y(j), 1-cgamma, z_asym, vareps);
end

%% Calculation
for i = 1:Nx
  for j = 1:MxG
    z = (exp(-(xmesh(i)-s*Y(j))^2) + exp(-(xmesh(i)+s*Y(j))^2));
    strain(i) = strain(i) + z*WT(j)*MF(j);
  end
end
strain = strain/2;
disp(['Current time: ', num2str(t), ' Negative Peak: ', num2str(min(min(strain)), '% 20.16f')])
toc

%% Plot
figure 
box on 
plot(xmesh, strain, '-o', 'LineWidth', 1, 'Markersize', 6);
xlabel('$$x$$', 'fontsize', 20, 'Interpreter', 'latex')
ylabel('strain', 'fontsize', 20, 'Interpreter', 'latex')
axis([-5 5 0 1])
drawnow;

% %% Save
% save varepsilon_cgamma_t varepsilon xmesh 
end