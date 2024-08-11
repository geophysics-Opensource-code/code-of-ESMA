function [err_inf, err_l2] = calculate_relative_error(sol_num, sol_ref, dx)

err_inf = max(abs(sol_num - sol_ref)/max(abs(sol_ref)));

err_l2  = sqrt(sum((sol_num - sol_ref).^2*dx));
err_l2 = err_l2/sqrt(sum((sol_ref).^2*dx));

% err_inf = max(abs(sol_num - sol_ref));
% err_l2  = sqrt(sum((sol_num - sol_ref).^2*dx));

end

