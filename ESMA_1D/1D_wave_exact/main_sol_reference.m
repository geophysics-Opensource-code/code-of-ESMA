%% Parameters
xmin     = -15;
xmax     = 15;
Nx       = 128;
tfin     = 8;
C        = 1;
%%%%%%%%%%%%%%%%
cgamma   = 0.1;
x_asym   = 1.852;
%%%%%%%%%%%%%%%%
% cgamma   = 0.05;
% x_asym   = 1.41;
%%%%%%%%%%%%%
% cgamma   = 0.01;
% x_asym   = 1.087;
%x_asym1   = sqrt(C)*tfin^(1-cgamma);
%%%%%%%%%%%%%
vareps = 10^(-7);
iflag = 0;

%% Grid mesh and initialization
dx     = (xmax - xmin)/Nx;
xmesh  = zeros(1, Nx);

% --- staggered-grid ---
for i = 1:Nx
  xmesh(i) = xmin + (i-1)*dx + dx/2;
end      

%% Calculation
if(tfin < 1)
  sol_ref = calculate_diffusive_wave_scaling(xmesh, tfin, C, cgamma, MxG, x_asym, vareps);

elseif(tfin >= 1)
  tic
  %--- Calculate the nodes and weights for the Gauss-Hermite quadrature ---
  if(iflag == 0)
    MxH = 200;
    sol_ref = calculate_diffusive_wave_Hermite(xmesh, tfin, C, cgamma, MxH, x_asym, vareps);
  
  %--- Calculate the nodes and weights for the Gauss-Legendre quadrature ---
  elseif(iflag == 1)  
    %% Preparing the Mainardi function
    MxG = 700;
    [Y, WT] = Laguerre_Gauss_FSU( MxG, 0, 1);
    MF = zeros(1, MxG);
    for j = 1:MxG
      MF(j) = calculate_Mainardi_function(Y(j), 1-cgamma, x_asym, vareps);
    end
  
    sol_ref = calculate_diffusive_wave_fixed(xmesh, tfin, C, cgamma, MF, Y, WT);
  
  elseif(iflag == 2)
    MxG = 800;
    s = sqrt(C)*tfin^(1-cgamma);
    [Y, WT] = Laguerre_Gauss_FSU( MxG, 0, s); 
    MF = zeros(1, MxG);
    for j = 1:MxG
      MF(j) = calculate_Mainardi_function(Y(j)/s, 1-cgamma, x_asym, vareps);
    end
  
    sol_ref = zeros(1, Nx);
    for i = 1:Nx
      for j = 1:MxG
        z = (exp(-(xmesh(i)-Y(j))^2) + exp(-(xmesh(i)+Y(j))^2));
        sol_ref(i) = sol_ref(i) + z*WT(j)*MF(j);
      end
    end
    sol_ref = sol_ref/(2*s);
    disp(['Current time: ', num2str(tfin), ' Positive Peak: ', num2str(max(sol_ref), '% 20.16f')])
  end
  toc
end

%% Plot
% figure 
% box on 
hold on
plot(xmesh, sol_ref, '-', 'LineWidth', 1, 'Markersize', 6);
xlabel('$$x$$', 'fontsize', 20, 'Interpreter', 'latex')
ylabel('Stress', 'fontsize', 20, 'Interpreter', 'latex')


