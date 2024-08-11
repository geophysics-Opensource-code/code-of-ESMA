function v_post = update_velocity_2D(v_pre, dsigma, dt,rou)

% --- Update the velocity (v) in leapfrog scheme

v_post = v_pre + dt* (1./rou) .*dsigma;

% v_post = v_pre + dt.*dsigma;

end

