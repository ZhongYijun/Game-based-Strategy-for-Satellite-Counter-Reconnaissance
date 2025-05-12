function cur_loss = Loss_cal(config, trans_mat, cur_state, cur_r, cur_ctl, J2L_mat, ...
    theta_max, unit_ball_sample, non_zero_count)
    % cur_non_zero_count = zeros(size(non_zero_count));
    % for k=1:config.num_rec_sat
    %     if norm(cur_state(k*config.orb_dim +(1:3),1))<=config.eff_obsv_distance
    %         cur_non_zero_count = cur_non_zero_count | (unit_ball_sample*cur_r(:,k)...
    %         >= cos(min(theta_max(k),0.5*pi - config.camera_view)));
    %     end
    % end
    % cur_non_zero_count = cur_non_zero_count & (unit_ball_sample*J2L_mat*config.a_sun...
    %                     <= cos(config.ang_sun));
    % non_zero_count = non_zero_count | cur_non_zero_count;
    % co_obser_area = 4*pi/config.max_num_ball_sample*sum(non_zero_count);
    V_R1 = 0;
    [fuel_loss_T, fuel_loss_R] = fuel_loss_cal(cur_ctl(1:config.ctl_dim,1),...
        cur_ctl(config.ctl_dim+1:2*config.ctl_dim,1), cur_ctl(2*...
        config.ctl_dim+1:end,1));
    V_T2 = config.wgt_rel_att_fuel*fuel_loss_T(1)^2 + fuel_loss_T(2)^2;
    V_R2 = fuel_loss_R^2;
    cur_payload_loss = payload_loss_cal(config, trans_mat, cur_r, J2L_mat, theta_max);
    
    V_T1 = config.wgt_rel_att_loss*cur_payload_loss;
    [V_T3, ~] = mission_loss(config, trans_mat, cur_r, J2L_mat, theta_max);
    cur_loss = zeros([2,1]);
    cur_loss(1) = config.wgt_T1*V_T1 + config.wgt_T2*V_T2 - config.wgt_T3*V_T3;
    cur_loss(2) = config.wgt_R1*V_R1 + config.wgt_R2*V_R2;

end