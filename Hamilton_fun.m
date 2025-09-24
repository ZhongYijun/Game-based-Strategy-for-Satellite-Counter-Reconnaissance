function [hamilton_] = Hamilton_fun(coeff, cur_state, cur_ctl, nxt_state,...
    J2L_mat, theta_max, config, type_value_fun, unit_ball_sample, non_zero_count)

[trans_mat, cur_r] = att_trans_mat(config, cur_state);
% theta_max = zeros([config.num_rec_sat, 1]);
% for j=1:config.num_rec_sat
%     theta_max(j) = acos(config.radius_tgt_sat/norm(cur_state((j-1)*...
%           config.orb_dim+config.att_dim+1:(j-1)*config.orb_dim+...
%           config.att_dim+3, 1)));
% end
hamilton_f = coeff'*(config.eta*basis_fun(config.wgt_mat*nxt_state) - ...
    basis_fun(config.wgt_mat*cur_state)) + Loss_cal(config, trans_mat,...
    cur_state, cur_r, cur_ctl, J2L_mat, theta_max, unit_ball_sample, non_zero_count);
if type_value_fun == "V_T"
    hamilton_ = hamilton_f(1);
else
    hamilton_ = hamilton_f(2);
end