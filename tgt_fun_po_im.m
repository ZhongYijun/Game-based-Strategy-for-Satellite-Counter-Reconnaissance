function [tgt_fun_value] = tgt_fun_po_im(coeff, cur_state, config, J2L_mat,...
    theta_max, cur_ctl, lst_ctl, type_value_fun, unit_ball_sample, non_zero_count)
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

nxt_state= Simulator(config, cur_state, cur_ctl);
loss_value = Loss_cal(config, trans_mat, cur_state, cur_r, cur_ctl, J2L_mat,...
    theta_max, unit_ball_sample, non_zero_count);

if type_value_fun == "V_T"
    tgt_fun_value = config.eta*coeff'*basis_fun(config.wgt_mat*nxt_state)...
        + loss_value(1);
else
    tgt_fun_value = config.eta*coeff'*basis_fun(config.wgt_mat*nxt_state)...
        + loss_value(2);
end

%% --- START OF SeBIL-GNE REGULARIZATION FIX ---
% Add the regularization penalty term: τ * ||u_current - u_anchor||^2
% This term "pulls" the new policy towards the last stable policy,
% preventing divergence.

if type_value_fun == "V_T"
    % Regularize only the Target's control
    player_ctl_current = cur_ctl(1 : 2*config.ctl_dim);
    player_ctl_anchor = lst_ctl(1 : 2*config.ctl_dim);
else % "V_R"
    % Regularize only the Recon's control
    player_ctl_current = cur_ctl(2*config.ctl_dim+1 : end);
    player_ctl_anchor = lst_ctl(2*config.ctl_dim+1 : end);
end

divergence_penalty = config.tau_regularization * norm(player_ctl_current - player_ctl_anchor);

%% Combine them into the final objective value
tgt_fun_value = tgt_fun_value + divergence_penalty;
% --- END OF FIX ---

end