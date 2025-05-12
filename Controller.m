function [ini_admissible_ctl] = Controller(config, cur_state,...
    trans_mat, cur_r, his_obs, cur_ctl, J2L_mat, unit_ball_sample, ...
    non_zero_count)
% Controller 
%% 首先来对于状态进行采样，采样方法为sobol序列采样
p1 = sobolset(3,'Leap',2);
p2 = sobolset(3,'Leap',3);
p3 = sobolset(3,'Leap',5);
p4 = sobolset(3,'Leap',7);

state_sample_stack = zeros([size(cur_state,1),config.num_sample]);
state_sample_stack(1:config.att_dim,:) = cur_state(1:config.att_dim) +...
    [0.1*p1(1:config.num_sample, :)';0.02*p2(1:config.num_sample, :)']; 

for i=1:config.num_rec_sat
    state_sample_stack(config.att_dim+((i-1)*config.orb_dim+1:i*config.orb_dim)...
        ,:) = cur_state(config.att_dim+((i-1)*config.orb_dim+1:i*config.orb_dim))...
        + [10*p3(1:config.num_sample, :)';0.1*p4(1:config.num_sample, :)'];
end
%% 接下来计算初始可行策略
[ini_admissible_ctl] = mpc_controller(config, cur_state, ...
    state_sample_stack, trans_mat, cur_r, his_obs, cur_ctl, J2L_mat, ...
    unit_ball_sample, non_zero_count);
%% 对于初始可行策略进行提升
% cur_ctl = qcqp_pi_controller(config, cur_state, state_sample_stack, his_obs,...
%     ini_admissible_ctl, J2L_mat);
% cur_ctl = ini_admissible_ctl(:,1);

end