function [nxt_costate] = costate_equation(config, cur_costate, cur_state,...
    his_obs, cur_ctl, theta_max)
%
nxt_state_diff_cur = jacobian(@(cur_state) Simulator(config, cur_state,...
        cur_ctl), cur_state);
loss_diff_state = gradient(@(cur_state) Loss_cal(config, cur_state,...
    cur_ctl, his_obs, theta_max), cur_state);
nxt_costate = nxt_state_diff_cur'\(cur_costate - loss_diff_state);

end