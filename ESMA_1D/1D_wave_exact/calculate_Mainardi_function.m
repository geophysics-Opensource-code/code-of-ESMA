function y = calculate_Mainardi_function(z, nu, z_asym, vareps)

%% calculate_Mainardi_function: calculate the Mainardi function adaptively through the Gauss-Laguerre rule.
%

%  Input:
%
%    z: number of points in the rule;
%
%    real nu, the exponent;
%

%% zero point
if( z == 0)
  y = 1/gamma(1-nu);
  
%% Laguerre-Gauss   
elseif(abs(z) <= z_asym)
  Mx = 100; 
  y0 = calculate_Wright_function(z, -nu, 1-nu, Mx, z_asym);
  y1 = -1;
  
  while(abs(y1 - y0) > vareps)
    y1 = y0;
    Mx = Mx + 50;
    y0 = calculate_Wright_function(z, -nu, 1-nu, Mx, z_asym);
  
    if(Mx > 3000)
      disp(['Fail to convergence in calculating the Mainardi function:  ', num2str(z), ' ', num2str(nu)])
      break;
    end
  end
  
  if(Mx < 3000)
    y = y0;
  else
    y = calculate_Wright_function(z, -nu, 1-nu, 1000, z_asym);
  end

%% Asymptotic expansion
else
  y = calculate_Mainardi_asymtotic(abs(z), nu);
end
  
end
  


