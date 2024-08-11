 %--------------ESMA2D-----------------------------
 %This program simulates the propagation of viscoacoustic waves in a 2D heterogeneous  model using the ESMA algorithm.
 % With the same velocity but different quality factors, different attenuation levels as well as weakly reflected waves can be observed
% Author: [Feng Wang]
%Require: data for wave velocity and the Q model to the 2D gas chimney model.
% Date: August 9, 2024
 clear
 close  all
 tic

xmin   = 0;
xmax   = 2000;
zmin   = 0;
zmax   = 2000;
LX     = xmax - xmin;
LZ     = zmax - zmin; 
Nx     = 200;
Nz     = 200;
My     = 8;

Q1     = 32;
Q2     = 32;
Q3     = 5;
Q4     = 5;
cgamma1 = (1/pi)*atan(1/Q1);%fractional order
cgamma2 = (1/pi)*atan(1/Q2);
cgamma3 = (1/pi)*atan(1/Q3);
cgamma4 = (1/pi)*atan(1/Q4);


dt     = 0.25e-3;
T      = 0.35;
Nstep  = fix(T/dt);

beta   = 0.05346*T + 0.209;%Scale factor
rou    = 2200*ones(Nx,Nz);

C1    = 2164.*ones(  Nx/2, Nz/2 );%
C2    = 2164.*ones( Nx/2, Nz/2 );
C3    = 2164.*ones( Nx/2, Nz/2 );
C4    = 2164.*ones( Nx/2, Nz/2 );

%------------Riker wavelet---------
fp = 25;
dr = 20*dt;

fp_REF = 1500;
W0   = 2*pi*fp;

C1    = C1.^2.*(cos(pi*cgamma1./2)).^2*W0.^(-2*cgamma1);
C2    = C2.^2.*(cos(pi*cgamma2./2)).^2*W0.^(-2*cgamma2);
C3    = C3.^2.*(cos(pi*cgamma3./2)).^2*W0.^(-2*cgamma3);
C4    = C4.^2.*(cos(pi*cgamma4./2)).^2*W0.^(-2*cgamma4);

tic
% --- Generation ---------
[node1, omega1] = Laguerre_Gauss_FSU(My, 4*cgamma1-1, beta);
delete('temp*')
Y1  = node1;
WT1 = omega1;

[node2, omega2] = Laguerre_Gauss_FSU(My, 4*cgamma2-1, beta);
delete('temp*')
Y2  = node2;
WT2 = omega2;


[node3, omega3] = Laguerre_Gauss_FSU(My, 4*cgamma3-1, beta);
delete('temp*')
Y3  = node3;
WT3 = omega3;

[node4, omega4] = Laguerre_Gauss_FSU(My, 4*cgamma4-1, beta);
delete('temp*')
Y4  = node4;
WT4 = omega4;

coeff1   = ones(  Nx/2, Nz/2 );
coeff2   = ones( Nx/2, Nz/2 );
coeff3   = ones(  Nx/2, Nz/2 );
coeff4   = ones( Nx/2, Nz/2 );


Y = zeros(My,Nx,Nz);
WT = zeros(My,Nx,Nz);

for ii = 1:My

Y (ii,:,:) = [Y1(ii)*coeff1,Y2(ii)*coeff2;Y3(ii)*coeff3,Y4(ii)*coeff4];

WT (ii,:,:) = [WT1(ii)*coeff1,WT2(ii)*coeff2;WT3(ii)*coeff3,WT4(ii)*coeff4];

end


% --- initialization -----------
dx      = (xmax - xmin)/Nx;
x1      =  linspace(xmin, xmax, Nx);
[xmesh,zmesh] = meshgrid(x1,(x1));
xmesh_shift = xmesh + dx/2;
zmesh_shift = zmesh + dx/2;

v1    =  zeros( Nx,Nx);
v3    =  zeros( Nx,Nx); 
phi   = zeros(My,Nx,Nx);
sigma = zeros(Nx,Nz);


s  = wavelet(fp,dt,T);

nn   = 6;
dxd  =  FDcoeffDx(nn);
ddz0 =  dxd';
ddz1 =  [dxd 0]';
ddx0 =  ddz0';
ddx1 =  ddz1';

pic_num = 1;   



for k = 1:Nstep-1 

%------------------------Stress excitation source loading------------------------ 

 F1 = 0;
 F2 = 0;

%----------------velocity--------------

  dsigmax = imfilter(sigma,ddx0)./dx;
  dsigmaz = imfilter(sigma,ddz0)./dx;

    v1 = update_velocity_2D(v1, dsigmax, dt,rou);
    v3 = update_velocity_2D(v3, dsigmaz, dt,rou);
    
%-----------------stress--------------

   dv1 = imfilter(v1,ddx1)./dx;
   dv3 =  imfilter(v3,ddz1)./dx;

    [sigma, phi] = update_stress_response_2D_changeQ(phi, Y,WT, dv1,dv3, xmesh,dt, C1,C2,C3,C4,...
        cgamma1,cgamma2,cgamma3,cgamma4,rou,F1,F2);
  
 sigma(3*Nz/8,Nx/2)  =  sigma(3*Nz/8,Nx/2) + s(k);


end
toc

figure(100)
pcolor(flipud(sigma))
shading interp;
colormap("bone");
colorbar;










