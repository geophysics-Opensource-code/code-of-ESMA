function vec = calculate_diffusive_wave_fixed(xmesh, t, C, cgamma, MF, Y, WT)
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
vec = zeros(1, Nx);

%% Calculation
s = sqrt(C)*t^(1-cgamma);
MxG = size(MF, 2);
for i = 1:Nx
  for j = 1:MxG
    z = (exp(-(xmesh(i)-s*Y(j))^2) + exp(-(xmesh(i)+s*Y(j))^2));
    vec(i) = vec(i) + z*WT(j)*MF(j);
  end
end
vec = vec/2;

for i = 1:Nx
  for j = 1:MxG
    z = (exp(-(xmesh(i)-s*Y(j))^2) + exp(-(xmesh(i)+s*Y(j))^2));
    vec(i) = vec(i) + z*WT(j)*MF(j);
  end
end
vec = vec/2;

disp(['Current time: ', num2str(t), ' Positive Peak: ', num2str(max(max(vec)), '% 20.16f')])

%% Plot
% figure 
% box on 
% plot(xmesh, strain, '-o', 'LineWidth', 1, 'Markersize', 6);
% xlabel('$$x$$', 'fontsize', 20, 'Interpreter', 'latex')
% ylabel('strain', 'fontsize', 20, 'Interpreter', 'latex')
% axis([xmin, xmax, 0, 1.1])
% drawnow;

% %% Save
% save varepsilon_cgamma_t varepsilon xmesh 
end