function dsigma = assemble_shifted_stress_derivate( sigma, L, M )

% --- Constant ---
b0 = 1i*pi/L;
b1 = 1i*pi/M;

% --- Fast Fourier transform ---
b = fft(sigma);

% --- Update mode ---
b(1:M/2+1) = 2*b0*b(1:M/2+1).*(0:1:M/2).*exp(b1*(0:1:M/2));
b(M/2+2:M) = 2*b0*b(M/2+2:M).*(-M/2+1:1:-1).*exp(b1*(-M/2+1:1:-1));

% --- Inverse Fourier transform
dsigma = ifft(b);
dsigma = real(dsigma);

end

