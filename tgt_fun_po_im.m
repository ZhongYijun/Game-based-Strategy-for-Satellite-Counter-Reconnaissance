function [tgt_fun_value] = tgt_fun_po_im(coeff, cur_state, config, his_obs,...
    J2L_mat, theta_max, cur_ctl, type_value_fun)
% 
[trans_mat, cur_r] = att_trans_mat(config, cur_state);
% if type_opt == "unconstraint"
% nxt_state_diff_cur = jocabian(@(cur_state) Simulator(config, cur_state,...
%     cur_ctl), cur_state);
% basis_fun_diff_state = jocabian(@(cur_state) basis_fun(cur_state), cur_state);
% loss_diff_state = gradient(@(cur_state) Loss_cal(config, trans_mat, cur_r,...
%     zeros(config.ctl_dim*(2+config.num_rec_sat), 1), his_obs), cur_state);
% 
% costate_diff = -basis_fun_diff_state'*coeff + loss_diff_state + config.eta*...
%     nxt_state_diff_cur'*cur_ctl(1);
% tgt_fun_value = costate_diff'*costate_diff;
% 
% else
nxt_state= Simulator(config, cur_state, cur_ctl);
cur_r = Loss_cal(config, trans_mat, cur_r, cur_ctl, his_obs, J2L_mat, theta_max);
if type_value_fun == "V_T"
    tgt_fun_value = config.eta*coeff'*basis_fun(config.wgt_mat*nxt_state)...
        + cur_r(1);
else
    tgt_fun_value = config.eta*coeff'*basis_fun(config.wgt_mat*nxt_state)...
        + cur_r(2);
end

% end
end