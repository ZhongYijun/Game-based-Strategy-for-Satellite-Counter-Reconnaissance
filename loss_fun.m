function [trun_v] = loss_fun(config, state, ctl_seq, his_obs, J2L_mat,...
    type_value_fun, unit_ball_sample, non_zero_count)
%
cur_state = state;
% cur_costate = costate;
V_T = 0;
V_R = 0;
% opts = optimoptions("fmincon", "Algorithm","sqp");

% cur_ctl_orb_rec = zeros([config.ctl_dim*config.num_rec_sat,1]);
for i = 1:config.T_sampling
    [trans_mat, cur_r] = att_trans_mat(config, cur_state);
    % his_obs = his_obs_update(config, cur_r, his_obs);
    % cur_ctl_att_tgt = ctl_seq([1:config.ctl_dim, i]);
    % cur_ctl_orb_tgt = ctl_seq([config.ctl_dim+1:2*config.ctl_dim, i]);

    % for j=1:config.num_rec_sat
    %     cur_ctl_orb_rec(:,j) = ctl_seq([(1+j)*config.ctl_dim+1: j*config.ctl_dim...
    %         +1, i]);
    % end
    cur_ctl = ctl_seq(:, i);
    theta_max = zeros([config.num_rec_sat, 1]);
    for j=1:config.num_rec_sat
        theta_max(j,1) = acos(config.radius_tgt_sat/norm(cur_state(...
            config.att_dim+(j-1)*config.orb_dim+(1:3), 1)));
    end
    
    ins_loss = Loss_cal(config, trans_mat, cur_state, cur_r, cur_ctl, ...
        J2L_mat, theta_max, unit_ball_sample, non_zero_count);
    V_T = V_T + ins_loss(1);
    V_R = V_R + ins_loss(2);
    for k=1:config.num_rec_sat
    if norm(cur_state(k*config.orb_dim +(1:3),1))>config.eff_obsv_distance
        V_R = V_R+config.wgt_dis_punish*norm(cur_state(k*config.orb_dim +(1:3),1))^2;
    elseif norm(cur_state(k*config.orb_dim +(1:3),1))<config.safe_distance
        V_R = V_R-config.wgt_dis_punish*norm(cur_state(k*config.orb_dim +(1:3),1))^2;
    end
    end
    cur_state = Simulator(config, cur_state, cur_ctl);
    % cur_costate = costate_equation(config, cur_costate, cur_state,...
    % his_obs, cur_ctl, theta_max);
end
cur_non_zero_count = zeros(size(non_zero_count));
for k=1:config.num_rec_sat
    if norm(cur_state(k*config.orb_dim +(1:3),1))<=config.eff_obsv_distance
        cur_non_zero_count = cur_non_zero_count | (unit_ball_sample*cur_r(:,k)...
        >= cos(min(theta_max(k),0.5*pi - config.camera_view)));
    end
end
cur_non_zero_count = cur_non_zero_count & (unit_ball_sample*J2L_mat*config.a_sun...
                        <= cos(config.ang_sun));
non_zero_count = non_zero_count | cur_non_zero_count;
co_obser_area = 4*pi/config.max_num_ball_sample*sum(non_zero_count);
V_T = V_T + config.wgt_T1_obs*co_obser_area;
V_R = V_R - config.wgt_R1_obs*co_obser_area;
if type_value_fun == "V_T"
    trun_v = V_T;
else
    trun_v = V_R;
end