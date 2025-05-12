function [co_obser_area] = tgt_seq_gen(config, cur_state, J2L_mat, theta_max,...
     tgt_r, type_value_fun, unit_ball_sample, non_zero_count_stack)
%
% 输入：当前状态: cur_state 历史观测区域近似: his_obs
% co_obser_area = integral2(@(psi,theta) tgt_state_generator(config, psi, theta,...
%     his_obs, cur_r, J2L_mat, theta_max, tgt_r), 0, pi,0, 2*pi,'AbsTol', 1e-6, ...
%     'RelTol', 1e-6);
if size(tgt_r,1) ~= 3
    tgt_r = tgt_r';
end
cur_non_zero_count = zeros(size(non_zero_count_stack));
for k=1:config.num_rec_sat
    if norm(cur_state(k*config.orb_dim +(1:3),1))<=config.eff_obsv_distance
       cur_non_zero_count = cur_non_zero_count | (unit_ball_sample*tgt_r(:,k)...
            >= cos(min(theta_max(k),0.5*pi - config.camera_view)));
    end
end
cur_non_zero_count = cur_non_zero_count & (unit_ball_sample*J2L_mat*config.a_sun...
                        <= cos(config.ang_sun));
non_zero_count = non_zero_count | cur_non_zero_count;
co_obser_area = 4*pi/config.max_num_ball_sample*sum(non_zero_count);


if type_value_fun == "V_R"
    co_obser_area = -co_obser_area;
end

end