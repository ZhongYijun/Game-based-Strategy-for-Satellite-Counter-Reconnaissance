function [non_zero_count, non_zero_count_stack] = non_zero_count_updater(config, ...
    theta_max, J2L_mat, unit_ball_sample, non_zero_count, non_zero_count_stack, ...
    cur_state, cur_r)
cur_non_zero_count = zeros(size(non_zero_count, 1), config.num_rec_sat);
for k=1:config.num_rec_sat
if norm(cur_state(k*config.orb_dim +(1:3),1)) <= config.eff_obsv_distance

    cur_non_zero_count(:,k) = cur_non_zero_count(:,k) | (unit_ball_sample*cur_r(:,k)...
                        >= cos(min(theta_max(k),0.5*pi - config.camera_view)));
    cur_non_zero_count(:,k) = cur_non_zero_count(:,k) & (unit_ball_sample*J2L_mat ...
                        *config.a_sun <= cos(config.ang_sun));
    non_zero_count_stack(:,k) = non_zero_count_stack(:,k) | cur_non_zero_count(:,k);
    non_zero_count = non_zero_count | cur_non_zero_count(:,k);


end
end
end