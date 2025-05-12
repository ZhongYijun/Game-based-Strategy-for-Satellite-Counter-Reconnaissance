function [V_T3, rel_ang] = mission_loss(config, trans_mat, cur_r, J2L_mat,...
    theta_max)
%% 函数指标用于反映目标卫星的任务连续性
% 输入：LVLH 轨道坐标系 to 体坐标系  trans_mat
% 输出：指标值V_T4, 当前摄像机指向与对地指向的误差角 rel_ang
rel_ang = config.q(:,1)'*trans_mat*config.a_mission;
cand_q_r =zeros([length(config.q),config.num_rec_sat]);
for i = 1:config.num_payload
    for j=1:config.num_rec_sat
        cand_q_r(i,j) = h1(config.q(:, i)'*trans_mat*cur_r(:,j), cos(min(...
            theta_max(j), (0.5*pi-config.camera_view))+config.safe_ang));
    end
end
V_T3 = (1-max(max(cand_q_r)))*h1(config.q(:,1)'*trans_mat*J2L_mat*config.a_sun,...
    -cos(config.ang_sun))*config.q(:,1)'*trans_mat*config.a_mission;

end