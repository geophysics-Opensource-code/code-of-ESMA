function y = calculate_Mainardi_asymtotic(z, nu)

%% calculate_Mainardi_asymtotic: calculate the first term in the asymptotic expansion of the Mainardi function

tz = z*nu;
a = 1/sqrt(2*pi*(1-nu));
b = (1 - nu)/nu;

y = a*tz.^((nu - 1/2)/(1 - nu)).*exp(-b*tz.^(1/(1-nu)));


end

