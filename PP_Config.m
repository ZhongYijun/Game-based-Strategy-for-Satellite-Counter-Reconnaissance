classdef PP_Config
    properties
        %% 系统动力学相关
        att_dim = 6                     % 姿态动力学维度 
        orb_dim = 6                     % 轨道动力学维度                                        
        att_ctl_dim = 3                 % 姿态发动机自由度
        orb_ctl_dim = 3                 % 轨道发动机自由度
        ctl_dim = 3
        mu_ = 3.986*10^14               % 地球引力常数
        MI = [611, -21, -21; -21, 607, -2; -21, -2, 666]
        % 目标卫星转动惯量
        m = 1185                        % 航天器质量
        orb_ele = [6778137; 0; 98*pi/180; pi/6; pi/2; 0]   
        % 虚拟航天器轨道六要素(近geo轨道)
        orb_w = 0;                      % 虚拟航天器轨道角速度
        J2 = 1.08*10^(-3);              % J2摄动量
        A_1 = zeros([3,3]);
        A_2 = zeros([3,3]);                          
        A_J2 = zeros([3,3]);
        B = [zeros([3,3]); eye(3)];     % 考虑J2摄动下的线性化轨道动力学方程                                                
        A_ZOH = zeros([6,6]);
        B_ZOH = zeros([6,3]);           % ZOH离散化后线性系统对应的系数矩阵
        %% 状态初值设置
        tgt_att_0 = [0;0;0;0;0;0]       % 目标卫星初始姿态状态 (6维)
        rel_orb_0 = zeros(6,1)          % 初始相对状态 (6维)
        r_e = 6.3781*10^6               % 地球轨道半径
        %% 控制器参数设置
        sample_time_span = 1            % 采样时间间隔
        control_time_span = 1           % 控制时间间隔
        %% MPC参数相关设置
        T_sampling = 5                  % MPC中中滚动时域步长
        T_control = 1                   % MPC中每步采取的控制量个数
        mpc_sample_time = 1             % MPC中采样时间间隔
        %% 博弈设置相关
        num_rec_sat = 1                 % 侦察卫星的个数
        num_payload = 3                 % 目标卫星载荷的个数
        tgt_att_u_max = 10              % 目标卫星姿态发动机最大力矩(nm)
        tgt_orb_u_max = 0.01            % 目标卫星轨道发动机最大推力加速度
        rec_orb_u_max = [0.01;0.01;0.01]     % 侦察卫星集群轨道发动机最大推力加速度
        q = [-1, -sqrt(3)/2, -sqrt(3)/2; 0, 1/2, -1/2; 0, 0, 0]
        % 定义相机矢量在本体坐标系中的位置矢量
        a_mission =[-1; 0; 0]           % 设置轨道坐标系下反映卫星任务连续性的指向
        % 此处以摄像头对地指向来反映卫星的任务连续性(轨道坐标系下)
        max_sim = 10000                 % 最大仿真次数
        camera_view = pi/60             % 设置相机观测角
        radius_tgt_sat = 6             % 设置目标卫星包络球半径
        a_sun = [1; 0; 0]              % 设置J2000坐标系下的太阳向量指向
        ang_sun = pi/6                  % 设置太阳规避角
        safe_distance = 2000            % 设置安全距离
        eff_obsv_distance = 10000       % 设置相机有效观测距离
        wgt_T1 = 100                    % 任务T1的相对权重
        wgt_T1_obs = 10000
        wgt_R1_obs = 10000
        wgt_T2 = 10                     % 任务T2的相对权重
        wgt_T3 = 10                     % 任务T3的相对权重
        wgt_T4 = 1                      % 任务T4的相对权重
        wgt_R1 = 100                    % 任务R1的相对权重
        wgt_R2 = 10                    % 任务R2的相对权重
        wgt_rel_att_fuel = 0.005;
        wgt_rel_att_loss = 1;
        %% ADP相关参数设置
        cst_flag = 1                    % 约束存在标识符（布尔变量）
        eta = 1-10^(-4)                 % 设置折扣因子
        num_sample = 0                  % 设置采样样本个数
        max_policy_iter = 20            % 设置最大的策略迭代步数
        err_tol = 0.01                  % 设置值函数终止迭代容许误差
        wgt_mat = eye(18)*diag([ones([6,1]);10^(-3)*ones([3,1]);...
            ones([3,1]);10^(-3)*ones([3,1]); ones([3,1])])
        max_trust_redius = 100          %设置信赖域半径上界
        wgt_dis_punish = 1e-4;
        %% 杂项
        max_num_ball_sample = 100000;
        safe_ang = pi/18;
        max_unchanged = 10;
    end
    methods 
        % 无参构造函数
        function obj = PP_Config()
           vel_tgt = sqrt(obj.mu_/obj.orb_ele(1));
           tgt_state = [obj.orb_ele(1); 0; 0; 0; vel_tgt; 0];
           rand_mat = rand(3,obj.num_rec_sat);
           r_vec = rand_mat ./ vecnorm(rand_mat,obj.num_rec_sat);
           rel_r = 3000*r_vec;
           rel_vel = zeros([3,obj.num_rec_sat]);
           for i=1:obj.num_rec_sat
               r_rec = vecnorm(tgt_state(1:3,1) + rel_r, 2);
               norm_rel_vel = sqrt(obj.mu_ ./ r_rec)-vel_tgt; 
               rel_vel(:, i) = norm_rel_vel(i)*cross(r_vec(:,i), [1;0;0]);
           end

           A = zeros([obj.orb_dim, obj.orb_dim]);
           A(1:3, 4:6) = eye(3);
           omega_ref = vel_tgt/obj.orb_ele(1);
           obj.orb_w = omega_ref;
           opts = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt',...
               'MaxFunctionEvaluations', 10000);
           obj.rel_orb_0(1:6, 1) = fsolve(@(cur_state) orb_ini(omega_ref,...
               cur_state), [0;0;-2500;0;0.5;0], opts);
           obj.rel_orb_0(1:6, 2) = fsolve(@(cur_state) orb_ini(omega_ref,...
                cur_state), [0;0;-2400;0;-0.5;0], opts);
           obj.rel_orb_0(1:6, 3) = fsolve(@(cur_state) orb_ini(omega_ref,...
                cur_state), [0;0;-2450;0;0;0], opts);
            % obj.rel_orb_0(1:6, 4) = fsolve(@(cur_state) orb_ini(omega_ref,...
           %     cur_state), [0;0;2400;0;0;0], opts);
           % obj.rel_orb_0(1:3,:) = [3000,3000;0,0;0,0];
           % obj.rel_orb_0(4:6,:) = [0,0.4;-0.4,0;6000*obj.orb_w,6000*obj.orb_w];
           A(4, 1) = 3*omega_ref^2; A(6, 6) = -omega_ref^2;
           A(4:6, 1:3) = A(4:6, 1:3)+ 6*obj.mu_*obj.J2*obj.r_e^2/...
               obj.orb_ele(1)^5*diag([1,-0.25,-0.75]);
           A(4:5, 4:5) = [0, 2*omega_ref; -2*omega_ref, 0];
           obj.A_ZOH = expm(A*obj.sample_time_span); 
           obj.B_ZOH = integral(@(t) expm(A.*t),0,obj.sample_time_span,...
               'ArrayValued',true)*obj.B;
        end
    end
end