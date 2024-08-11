function [node, omega] = Laguerre_Gauss_FSU(M, alpha, beta)
% Generate points for Laguerre-Gauss integration   产生拉盖尔-高斯积分的点

gen_laguerre_rule ( M, alpha, 0, 1, 'temp' )
node   = load('temp_x.txt');
weight = load('temp_w.txt');
delete('temp*')

omega = zeros(M, 1);
syms t
F = laguerreL(M+1, alpha, t);%LaguerreL函数
F = t.*exp(t)./(F.^2);
for i = 1:M
  t = node(i);
  omega(i) = vpa(subs(F));
end

% --- factorial ---
t = node(1);
F = laguerreL(M+1, alpha, t);
fac = weight(1)*(vpa(subs(F))^2)/node(1);%vpa计算符号的变量和函数的值%subs将符号表达式中的某位符号变量替换为指定的新变量

omega = omega*fac;
omega = double(omega);

%% final step: scaling%缩放
node = node/beta;
omega = omega/((beta)^(alpha+1));

end