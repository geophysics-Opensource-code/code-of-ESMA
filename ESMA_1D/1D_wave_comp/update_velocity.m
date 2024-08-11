function v_post = update_velocity(v_pre, dsigma, mesh_shift, dt)

% --- Update the velocity (v) in leapfrog scheme

v_post = v_pre + dt*1./mass_function(mesh_shift).*dsigma + dt*source_function(mesh_shift);

end

