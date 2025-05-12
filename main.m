%% MIAN entrance of the payload protection simulation
%% 参数初始化
clear;clc;
config = PP_Config();
state_seq_stack = zeros([config.att_dim + config.num_rec_sat*config.orb_dim,...
    config.max_sim]);
%% 下述创建的stack每列对应一个状态的损失
ctl_seq_stack = zeros([(2+config.num_rec_sat)*config.ctl_dim, config.max_sim]);
cur_ctl = zeros([(2+config.num_rec_sat)*config.ctl_dim, config.num_sample+1]);
payload_loss_stack = pi*ones([config.num_payload, config.max_sim]); 
co_vis_area_stack = zeros([1, config.max_sim]);
eny_consum_stack = zeros([2+config.num_rec_sat, config.max_sim]);
rel_dis = zeros([config.num_rec_sat, config.max_sim]); 
mission_loss_stack = zeros([1, config.max_sim]);
cur_loss_stack = zeros([2, config.max_sim]);
r_stack = zeros([config.num_rec_sat, config.max_sim]);
cur_state = [config.tgt_att_0; reshape(config.rel_orb_0, config.orb_dim...
    *config.num_rec_sat,1)];
[trans_mat, cur_r] = att_trans_mat(config, cur_state);
state_stack = [];
his_obs = zeros(4, config.num_rec_sat);
his_obs(2:end, :) = cur_r;
J2L_mat = zeros([3,3]);
unit_ball_sample = unit_ball_sample_gen(config);
non_zero_count = zeros([config.max_num_ball_sample, 1]);
non_zero_count_stack = zeros([config.max_num_ball_sample, config.num_rec_sat]);
r_rec = cur_r; r_tgt = cur_r;
%% 进行仿真
for i = 1:config.max_sim
    tic
    theta = config.orb_ele(6)+i*config.sample_time_span*config.orb_w; 
    % 更新J2000坐标系到LVLH坐标系的转换矩阵
    J2L_mat(1:3, :) = [cos(theta)*cos(config.orb_ele(4))-sin(theta)*sin(...
        config.orb_ele(4))*cos(config.orb_ele(3)), cos(theta)*sin(...
        config.orb_ele(4))+sin(theta)*cos(config.orb_ele(4))*cos( ...
        config.orb_ele(3)), sin(theta)*sin(config.orb_ele(3)); -sin(theta)...
        *cos(config.orb_ele(4))-cos(theta)*sin(config.orb_ele(4))*cos(...
        config.orb_ele(3)), -sin(theta)*sin(config.orb_ele(4))+cos(theta)*...
        cos(config.orb_ele(4))*cos(config.orb_ele(3)), cos(theta)*sin(...
        config.orb_ele(3)); sin(config.orb_ele(3))*sin(config.orb_ele(4)),...
        -sin(config.orb_ele(3))*cos(config.orb_ele(4)), cos(config.orb_ele(3))];
    theta_max = zeros([config.num_rec_sat, 1]);
    for j=1:config.num_rec_sat
            theta_max(j,1) = acos(config.radius_tgt_sat/norm(cur_state(...
           j*config.orb_dim+(1:config.ctl_dim), 1)));
    end
    
    % 设计控制器来保证性能
    [non_zero_count, non_zero_count_stack] = non_zero_count_updater(config, ...
        theta_max, J2L_mat, unit_ball_sample, non_zero_count, ...
        non_zero_count_stack, cur_state, cur_r);

    cur_ctl = Controller(config, cur_state, trans_mat, cur_r,...
         his_obs, cur_ctl, J2L_mat, unit_ball_sample, non_zero_count);
    if i>=config.max_unchanged && max(abs(mean(co_vis_area_stack(:,i+ ...
       (-config.max_unchanged+1:0)))-co_vis_area_stack(:,i+(-config.max_unchanged+1:0))))< 1e-4  
        cur_ctl(config.ctl_dim+(1:config.ctl_dim),1) = zeros([config.orb_dim, 1]);
    end

    % 显示当前状态并打印当前时间状态
    disp("current time:")
    disp(i)
    disp("current control vector:")
    disp(cur_ctl)
    disp("current relative distance vector:")
    disp(cur_r)
    
    for j = 1:config.num_payload
        for k=1:config.num_rec_sat
            if cur_r(:, k)'*J2L_mat*config.a_sun <= cos(config.ang_sun)
                payload_loss_stack(j, i) = min(payload_loss_stack(j, i), ...
                    acos(config.q(:,j)'*trans_mat*cur_r(:, k)));
            end
        end
    end

    co_vis_area_stack(:, i) = 4*pi/config.max_num_ball_sample*sum(...
                    non_zero_count);

    for j = 1:2+config.num_rec_sat
        eny_consum_stack(j, i) = norm(cur_ctl(config.ctl_dim*(j-1)+...
            (1:config.ctl_dim),1));
    end
    
    mission_loss_stack(:, i) = acos(config.q(:,1)'*trans_mat*config.a_mission);

    if  co_vis_area_stack(:, i) + 2*config.num_rec_sat*pi*(1-cos(...
            config.safe_ang)) < 4*pi
        
       cur_ctl(config.ctl_dim+(1:config.ctl_dim), 1) = zeros(config.ctl_dim, 1);
    end
    cur_state = Simulator(config, cur_state, cur_ctl(:,1));
    [trans_mat, cur_r] = att_trans_mat(config, cur_state);
    state_stack = [state_stack, cur_state];
    toc
end
set(0,'defaultfigurecolor','w')
figure(1)
plot(co_vis_area_stack, 'Linewidth', 2)
xlabel('Time(s)')
ylabel('Co-visible area')
title("Area of co-visible area of RS")
figure(2)
hold on 
for i=1:config.num_payload
    plot(payload_loss_stack(i, :),  'Linewidth', 1.5)
end
plot(ones([1,config.max_sim])*min(max(theta_max),0.5*pi - config.camera_view),...
    'Linewidth', 1)
hold off
xlabel('Time(s)')
ylabel('Payload loss(rad)')
title("Angle between payload and camera direction of RS")
legend(["payload 1", "payload 2", "payload 3"])
figure(3)
hold on
for i=2:2+config.num_rec_sat
    plot(eny_consum_stack(i, :), 'Linewidth', 1.5)
end
hold off
xlabel('Time(s)')
ylabel('Fuel loss')
title("Fuel loss of target satellite and reconnassiance satellite")
legend(["orb_TS", "orb_RS1", "orb_RS2", "orb_RS3"],'Interpreter','none')
figure(4)
plot(mission_loss_stack, 'Linewidth', 1.5)
xlabel('Time(s)')
ylabel('mission loss(rad)')
title("Angle between camera direction of TS and mission direction")
