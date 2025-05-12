function [payload_value] = payload_loss_cal(config, trans_mat, cur_r, J2L_mat, theta_max)
%% 函数用于计算载荷调整损失，当载荷位于安全位置外时，函数值不为0
payload_value = 0;
for j=1:config.num_payload
for i=1:config.num_rec_sat
    payload_value = payload_value + h1(config.q(:,j)'*trans_mat*cur_r(:,i),...
        cos(min(theta_max(i), (0.5*pi-config.camera_view)))).*h2(cur_r(:,i)'*...
        J2L_mat*config.a_sun, cos(config.ang_sun));
end
end

