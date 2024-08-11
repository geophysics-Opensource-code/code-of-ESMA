function  [sigma_post, phi_post] = update_stress_response_2D_Marmus (phi_pre, Y,wt, dv1,dv3, M,N, dt, C_gamma,rou,F1,F2,sigma_x,sigma_z)
% --- Update the stress (sigma) and response functional (phi) in leapfrog scheme

My = size(Y, 1);
sigma_post = zeros(N, M);
phi_post   = zeros(My, N,M);

% --- Constant ---

sigma_xz  = sigma_x+sigma_z;%PML衰减系数矩阵

  FS      = dt*(F1+F2); %外力做功

% --- Update the response functional ---
for i = 1:My
    
  TT  = exp( -squeeze(Y(i,:,:)).^2*dt).*squeeze( phi_pre(i,:,:) ) + ( 1./squeeze( Y(i,:,:) ).^2 ).*...
                    (1 - exp(-squeeze(Y(i,:,:)).^2*dt) ).*(dv1+dv3);

   phi_post(i, :,:) = TT - dt*sigma_xz.*TT;
  
end

% --- Update the stress ---



for i = 1:My
  sigma_post = sigma_post + squeeze(wt(i,:,:)).*squeeze( phi_post(i,:,:) );
end

sigma_post =( C_gamma.*rou ).*sigma_post + FS;

sigma_post  =   sigma_post - dt*sigma_xz.* sigma_post;

end

