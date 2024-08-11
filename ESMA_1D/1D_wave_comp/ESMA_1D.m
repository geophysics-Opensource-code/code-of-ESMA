%% plot_temporal_convergence: plot the convergence in operator splitting

% --- Coefficients ---
xmin   = -15;
xmax   =  15;
L      = xmax - xmin;
Nx     = 128;
dx     = L/Nx;

cgamma = 0.01;
    dt = 0.0001;
T      = 8;
C      = 1;
beta   = 1;

    % Laguerre-Gauss nodes and weights
    My = 32;

    [node, omega] = Laguerre_Gauss_FSU(My, 4*cgamma-1, beta);
    delete('temp*')
    Y  = node;
    WT = omega;
   

     %--- Initialization ---
    [mesh, mesh_shift, v, sigma, phi] = initialization(xmin, xmax, Nx, Y, cgamma);  
    Nstep  = fix(T/dt);

    tic
    % --- Time evolution ---
    for l = 1:Nstep
      % --- Half-step update of velocity --
      dsigma = assemble_shifted_stress_derivate( sigma, L, Nx );
      v      = update_velocity(v, dsigma, mesh_shift, dt/2);

      % --- Full-step update of stress
      dv = assemble_shifted_velocity_derivate( v, L, Nx );
      [sigma, phi] = update_stress_response(phi, Y, WT, dv, mesh, dt, C, cgamma);

      % --- Half-step update of velocity --
      dsigma = assemble_shifted_stress_derivate( sigma, L, Nx );
      v      = update_velocity(v, dsigma, mesh_shift, dt/2);
    end
    toc
    
    % --- calculate error ---
  
 plot(mesh_shift,v, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);









 

