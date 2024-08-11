function y = calculate_Wright_function(z, lambda, mu, Mx, z_asym)

%% calculate_Wright_function: calculate the Wright function through the Gauss-Laguerre rule.
%  Input:
%
%    z: number of points in the rule;
%
%    real nu, the exponent;
%
eps = 0.99;

%% case 1 : Gaussian
if(lambda == -1/2)
  gen_laguerre_rule ( Mx, -mu, 0, 1, 'genLG' )
  R  = load('genLG_x.txt');
  WT = load('genLG_w.txt');

  Y = func1(R, z, lambda, mu);
  phi = sum(WT.*Y);
  y   = phi/pi;
  delete('genLG*')
    
%% case 2: subdiffusion 
elseif(lambda > -1/2)
  if( z < 0 )
    y = calculate_Wright_function(-z, lambda, mu, Mx, z_asym);

  elseif(z > 0)
    s1 = (mu-1)/lambda - 1;
    s2 = cos(lambda*pi);
    gen_laguerre_rule ( Mx, s1, 0, s2, 'genLG' )
    R  = load('genLG_x.txt');
    WT = load('genLG_w.txt');

    Y = func2(R, z, lambda, mu);
    phi = sum(WT.*Y);
    y   = phi/pi;
    delete('genLG*')  
  end
  
%% case 3: superdiffusion  
elseif(lambda < -1/2 )
  if(z < 0)
    y = calculate_Wright_function(-z, lambda, mu, Mx, z_asym);

  elseif(z < (-1/cos(lambda*pi))*eps)
    s1 = -mu;
    s2 = (1 + z*cos(lambda*pi));
    
    gen_laguerre_rule ( Mx, s1, 0, s2, 'genLG' )
    R  = load('genLG_x.txt');
    WT = load('genLG_w.txt');

    Y = func3(R, z, lambda, mu);
    phi = sum(WT.*Y);
    y   = phi/pi;
    delete('genLG*') 

  elseif(z < -1/cos(lambda*pi))
    s1 = -mu;
    s2 = (1 + z*cos(lambda*pi));
    s2 = s2 + 0.05;
    
    gen_laguerre_rule ( Mx, s1, 0, s2, 'genLG' )
    R  = load('genLG_x.txt');
    WT = load('genLG_w.txt');

    Y = func5(R, z, lambda, mu, 0.05);
    phi = sum(WT.*Y);
    y   = phi/pi;
    delete('genLG*') 
    
  elseif(z >= -1/cos(lambda*pi))
    if(abs(z) > z_asym)
      Mx_asym = 100;
      Y = calculate_Wright_asymptotic(abs(z), lambda, mu, Mx_asym);
      y = sum(Y);
      
    else  
      if(lambda <= -0.9)
        Mx = 200; % lambda is small, Mx should not be very large for stability
      end
%       tz = -z*cos(lambda*pi);
%       y0 = tz^((1-lambda)/(1+lambda));
%       I1 = integral(@(r) func5(r, z, lambda, mu), 0, y0);
%       
      s1 = (mu-1)/lambda - 1;
      s2 = -z*cos(lambda*pi) - 1;
    
      gen_laguerre_rule ( Mx, s1, 0, s2, 'genLG' )
      R  = load('genLG_x.txt');
      WT = load('genLG_w.txt');

      Y = func4(R, z, lambda, mu, 0);
      I2 = sum(WT.*Y);
      y   = I2/pi;
      delete('genLG*') 
    end
  end
else
  y = 0;
end

end

%%
function y = func1(r, z, lambda, mu)

y = exp(-z.*r.^-lambda.*cos(lambda*pi)).*sin(-z.*r.^(-lambda).*sin(pi*lambda) + pi*mu);

end

%%
function y = func2(r, z, lambda, mu)

s1 = (mu - 1)/lambda - 1;
y = -(z^(-(s1+1)).*exp(-(z./r).^(1/lambda)).*sin(-r.*sin(pi*lambda) + pi*mu))./(lambda);

end

%%
function y = func3(r, z, lambda, mu)

y = exp(z.*(r-r.^-lambda).*cos(lambda*pi)).*sin(-z.*r.^(-lambda).*sin(pi*lambda) + pi*mu);

end

%%
function y = func4(r, z, lambda, mu, y0)

s1 = (mu - 1)/lambda - 1;
tz = z*(-cos(lambda*pi));
r = r + y0;
y = -(tz^(-(s1+1))*exp(-(1-tz)*y0).*exp(-(tz./r).^(1/lambda) + tz.*r).*sin(r.*tan(pi*lambda) + pi*mu))./(lambda);

end

%%
function y = func5(r, z, lambda, mu, eps)

y = exp(z.*(eps+(r-r.^-lambda).*cos(lambda*pi))).*sin(-z.*r.^(-lambda).*sin(pi*lambda) + pi*mu);

end




