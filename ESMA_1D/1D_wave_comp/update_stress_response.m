function [sigma_post, phi_post] = update_stress_response(phi_pre, Y, wt, dv, mesh, dt, C, cgamma)
% --- Update the stress (sigma) and response functional (phi) in leapfrog scheme
M = size(mesh, 2);
My = size(Y, 1);
sigma_post = zeros(1, M);
phi_post   = zeros(My, M);

% --- Constant ---
C_gamma = C*2*sin(2*pi*cgamma)/pi;

% --- Update the response functional ---
for i = 1:My
  phi_post(i, :) = exp(-Y(i)^2*dt)*phi_pre(i, :) + (1/Y(i)^2)*(1 - exp(-Y(i)^2*dt))*dv;
end

% --- Update the stress ---
for i = 1:My
  sigma_post = sigma_post + wt(i)*phi_post(i, :);
end
sigma_post = C_gamma*mass_function(mesh).*sigma_post;

end

