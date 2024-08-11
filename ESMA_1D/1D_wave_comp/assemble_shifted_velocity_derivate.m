function dv = assemble_shifted_velocity_derivate( v, L, M )

% --- Constant ---
a0 = 1i*pi/L;
a1 = 1i*pi/M;

% --- Fast Fourier transform ---
a = fft(v);

% --- Update mode ---
a(1:M/2+1) = 2*a0*a(1:M/2+1).*(0:1:M/2).*exp(-a1*(0:1:M/2));
a(M/2+2:M) = 2*a0*a(M/2+2:M).*(-M/2+1:1:-1).*exp(-a1*(-M/2+1:1:-1));

% --- Inverse Fourier transform
dv = ifft(a);
dv = real(dv);

end

