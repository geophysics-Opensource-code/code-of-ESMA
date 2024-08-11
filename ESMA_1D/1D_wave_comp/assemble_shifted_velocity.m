function vsh = assemble_shifted_velocity( v, L, M )

% --- Constant ---
a0 = 1i*pi/L;
a1 = 1i*pi/M;

% --- Fast Fourier transform ---
a = fft(v);

% --- Update mode ---
a(1:M/2+1) = a(1:M/2+1).*exp(-a1*(0:1:M/2));
a(M/2+2:M) = a(M/2+2:M).*exp(-a1*(-M/2+1:1:-1));

% --- Inverse Fourier transform
vsh = ifft(a);
vsh = real(vsh);

end

