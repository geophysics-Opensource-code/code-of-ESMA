function v_post = update_velocity_PML(v_pre, dsigma, dt,rou,sigma_xz)

% --- Update the velocity (v) in leapfrog scheme

v_post = v_pre + dt* (1./rou) .*dsigma;

v_post =   v_post  -dt*sigma_xz.*v_post;

% v_post = v_pre + dt.*dsigma;

end

