function [xmesh, xmesh_shift, v, sigma, phi] = initialization(xmin, xmax, Nx, Y, cgamma)
%% initialization: Initialization of grid mesh and all variables 

% --- grid mesh ---
dx   = (xmax - xmin)/Nx;
xmesh       = zeros(1, Nx);
xmesh_shift = zeros(1, Nx);
for i = 1:Nx
  xmesh(i)       = xmin + (i-1)*dx;
  xmesh_shift(i) = xmin + (i-1)*dx + dx/2;
end

% --- initialization ---
v = exp(-xmesh_shift.*xmesh_shift);
sigma = zeros(1, Nx);

My  = size(Y, 1);
phi = zeros(My, Nx);
for j = 1:My
  c = gamma(1-2*cgamma)*exp(-Y(j)^2);
  phi(j, :) = c*sigma;
end

end

