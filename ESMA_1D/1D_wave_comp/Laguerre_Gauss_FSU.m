function [node, omega] = Laguerre_Gauss_FSU(M, alpha, beta)
% Generate points for Laguerre-Gauss integration

gen_laguerre_rule ( M, alpha, 0, 1, 'temp' )
node   = load('temp_x.txt');
weight = load('temp_w.txt');
delete('temp*')

omega = zeros(M, 1);
syms t
F = laguerreL(M+1, alpha, t);
F = t.*exp(t)./(F.^2);
for i = 1:M
  t = node(i);
  omega(i) = vpa(subs(F));
end

% --- factorial ---
t = node(1);
F = laguerreL(M+1, alpha, t);
fac = weight(1)*(vpa(subs(F))^2)/node(1);

omega = omega*fac;
omega = double(omega);

%% final step: scaling
node = node/beta;
omega = omega/((beta)^(alpha+1));

end