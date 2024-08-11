function  [sigma_post, phi_post] = update_stress_response_2D_changeQ(phi_pre, Y,wt, dv1,dv3, xmesh,dt, C1,C2,C3,C4,...
    cgamma1,cgamma2,cgamma3,cgamma4,rou,F1,F2)
% --- Update the stress (sigma) and response functional (phi) in leapfrog scheme
[M ,N] = size(xmesh);
My = size(Y, 1);
sigma_post = zeros(M, N);
phi_post   = zeros(My, M,N);

% --- Constant ---

C_gamma1 = C1*2.*sin(2*pi*cgamma1)./pi;
C_gamma2 = C2*2.*sin(2*pi*cgamma2)./pi;
C_gamma3 = C3*2.*sin(2*pi*cgamma3)./pi;
C_gamma4 = C4*2.*sin(2*pi*cgamma4)./pi;

C_gamma = [C_gamma1,C_gamma2;C_gamma3,C_gamma4];


FS  = dt*(F1+F2);

% --- Update the response functional ---
for i = 1:My
    
  phi_post(i, :,:) = exp(-squeeze(Y(i,:,:)).^2*dt).*squeeze( phi_pre(i,:,:) ) + ( squeeze( 1./Y(i,:,:) ).^2 ).*...
      (1 - exp(-squeeze(Y(i,:,:)).^2*dt) ).*(dv1+dv3);
  
end

% --- Update the stress ---

for i = 1:My

  sigma_post = sigma_post + squeeze(wt(i,:,:)).*squeeze( phi_post(i,:,:) );
  
end

sigma_post =( C_gamma.*rou ).*sigma_post + 0.5*FS ;


end

