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
J2L_mat = zeros([3,3]);
unit_ball_sample = unit_ball_sample_gen(config);
non_zero_count = zeros([config.max_num_ball_sample, 1]);
r_rec = cur_r; r_tgt = cur_r;
%% 进行仿真
mean_time = 0;

for i = 1:config.max_sim

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
    non_zero_count = non_zero_count_updater(config, theta_max, J2L_mat, ...
        unit_ball_sample, non_zero_count, cur_state, cur_r);
    tic;
    cur_ctl = Controller(config, cur_state, cur_ctl, J2L_mat, unit_ball_sample,...
        non_zero_count);
    during_time = toc;
    if i>=config.max_unchanged && max(abs(mean(co_vis_area_stack(:,i+ ...
       (-config.max_unchanged+1:0)))-co_vis_area_stack(:,i+(-config.max_unchanged+1:0))))< 1e-4  
        cur_ctl(config.ctl_dim+(1:config.ctl_dim),1) = zeros([config.ctl_dim, 1]);
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
    mean_time = during_time/i + (i-1)*mean_time/i;
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

% function main()
% %% MIAN entrance of the payload protection simulation
% % =========================================================================
% %         最终修复版 - 最小改动以支持运行N次
% % =========================================================================
% % 说明：在之前版本基础上，用一个简单的for循环包裹了核心逻辑，
% % 以实现多次运行。每次运行的结果将保存到带序号的文件中。
% %==========================================================================
% 
% % 【新增改动】设置您希望运行的总次数
% num_runs = 3; 
% 
% % 【新增改动】创建一个cell数组来存储每次运行的结果（可选，但推荐）
% all_run_results = cell(1, num_runs);
% persistent last_cpu_usages;
% % 【新增改动】开始N次运行的大循环
% for run_index = 1:num_runs
% 
%     fprintf('\n\n=======================================================');
%     fprintf('\n          开始执行第 %d / %d 次仿真运行\n', run_index, num_runs);
%     fprintf('=======================================================\n\n');
% 
%     %% 参数初始化和变量定义 (此部分代码在每次循环开始时都会重新执行)
%     % clear; clc; % 在函数内部，clear和clc通常不是必需的，因为作用域是独立的
%     config = PP_Config(); 
% 
%     sim_start_time = datetime('now');
% 
%     current_cpu_usage = 0;
%     last_cpu_usages = [];
% 
%     old_timer = timerfind('Name', 'CPUMonitorTimer');
%     if ~isempty(old_timer)
%         stop(old_timer);
%         delete(old_timer);
%     end
% 
%     cpu_monitor = timer(...
%         'ExecutionMode', 'fixedRate', ...
%         'Period', 1, ...
%         'TimerFcn', @update_cpu_usage_robust, ...
%         'StartDelay', 0, ...
%         'Name', 'CPUMonitorTimer');
% 
%     performance_metrics = struct();
%     performance_metrics.sim_time = zeros(1, config.max_sim);
%     performance_metrics.memory_usage_MB = zeros(1, config.max_sim);
%     performance_metrics.cpu_usage = zeros(1, config.max_sim);
% 
%     %% 仿真状态变量初始化
%     state_seq_stack = zeros([config.att_dim + config.num_rec_sat*config.orb_dim, config.max_sim]);
%     ctl_seq_stack = zeros([(2+config.num_rec_sat)*config.ctl_dim, config.max_sim]);
%     cur_ctl = zeros([(2+config.num_rec_sat)*config.ctl_dim, config.num_sample+1]);
%     payload_loss_stack = pi*ones([config.num_payload, config.max_sim]);
%     co_vis_area_stack = zeros([1, config.max_sim]);
%     eny_consum_stack = zeros([2+config.num_rec_sat, config.max_sim]);
%     mission_loss_stack = zeros([1, config.max_sim]);
%     cur_state = [config.tgt_att_0; reshape(config.rel_orb_0, config.orb_dim*config.num_rec_sat,1)];
%     [trans_mat, cur_r] = att_trans_mat(config, cur_state);
%     state_stack = [];
%     his_obs = zeros(4, config.num_rec_sat);
%     his_obs(2:end, :) = cur_r;
%     J2L_mat = zeros([3,3]);
%     unit_ball_sample = unit_ball_sample_gen(config);
%     non_zero_count = zeros([config.max_num_ball_sample, 1]);
%     non_zero_count_stack = zeros([config.max_num_ball_sample, config.num_rec_sat]);
% 
%     %% 进行仿真
%     disp('仿真开始...');
%     start(cpu_monitor);
%     cleanupObj = onCleanup(@() stop_and_delete_timer(cpu_monitor));
% 
%     for i = 1:config.max_sim
%         iter_start_time = tic;
% 
%         % --- 核心仿真循环逻辑 (无改动) ---
%         theta = config.orb_ele(6)+i*config.sample_time_span*config.orb_w; 
%         J2L_mat(1:3, :) = [cos(theta)*cos(config.orb_ele(4))-sin(theta)*sin(config.orb_ele(4))*cos(config.orb_ele(3)), cos(theta)*sin(config.orb_ele(4))+sin(theta)*cos(config.orb_ele(4))*cos(config.orb_ele(3)), sin(theta)*sin(config.orb_ele(3)); -sin(theta)*cos(config.orb_ele(4))-cos(theta)*sin(config.orb_ele(4))*cos(config.orb_ele(3)), -sin(theta)*sin(config.orb_ele(4))+cos(theta)*cos(config.orb_ele(4))*cos(config.orb_ele(3)), cos(theta)*sin(config.orb_ele(3)); sin(config.orb_ele(3))*sin(config.orb_ele(4)), -sin(config.orb_ele(3))*cos(config.orb_ele(4)), cos(config.orb_ele(3))];
%         theta_max = zeros([config.num_rec_sat, 1]);
%         for j=1:config.num_rec_sat
%             theta_max(j,1) = acos(config.radius_tgt_sat/norm(cur_state(j*config.orb_dim+(1:config.ctl_dim), 1)));
%         end
%         [non_zero_count, non_zero_count_stack] = non_zero_count_updater(config, theta_max, J2L_mat, unit_ball_sample, non_zero_count, non_zero_count_stack, cur_state, cur_r);
%         cur_ctl = Controller(config, cur_state, trans_mat, cur_r, his_obs, cur_ctl, J2L_mat, unit_ball_sample, non_zero_count);
%         if i>=config.max_unchanged && max(abs(mean(co_vis_area_stack(:,i+(-config.max_unchanged+1:0)))-co_vis_area_stack(:,i+(-config.max_unchanged+1:0))))< 1e-4  
%             cur_ctl(config.ctl_dim+(1:config.ctl_dim),1) = zeros([config.ctl_dim, 1]);
%         end
%         if mod(i, 10) == 0, fprintf('当前迭代: %d / %d\n', i, config.max_sim); end
%         for j = 1:config.num_payload
%             for k=1:config.num_rec_sat
%                 if cur_r(:, k)'*J2L_mat*config.a_sun <= cos(config.ang_sun)
%                     payload_loss_stack(j, i) = min(payload_loss_stack(j, i), acos(config.q(:,j)'*trans_mat*cur_r(:, k)));
%                 end
%             end
%         end
%         co_vis_area_stack(:, i) = 4*pi/config.max_num_ball_sample*sum(non_zero_count);
%         for j = 1:2+config.num_rec_sat
%             eny_consum_stack(j, i) = norm(cur_ctl(config.ctl_dim*(j-1)+(1:config.ctl_dim),1));
%         end
%         mission_loss_stack(:, i) = acos(config.q(:,1)'*trans_mat*config.a_mission);
%         if  co_vis_area_stack(:, i) + 2*config.num_rec_sat*pi*(1-cos(config.safe_ang)) < 4*pi
%            cur_ctl(config.ctl_dim+(1:config.ctl_dim), 1) = zeros(config.ctl_dim, 1);
%         end
%         cur_state = Simulator(config, cur_state, cur_ctl(:,1));
%         [trans_mat, cur_r] = att_trans_mat(config, cur_state);
%         state_stack = [state_stack, cur_state];
% 
%         performance_metrics.sim_time(i) = toc(iter_start_time);
%         performance_metrics.memory_usage_MB(i) = get_process_resident_memory() / (1024^2);
%         performance_metrics.cpu_usage(i) = current_cpu_usage;
%     end
% 
%     disp('仿真结束.');
% 
%     %% 性能评估与结果保存
%     disp('正在计算性能指标并保存结果...');
%     sim_end_time = datetime('now');
%     performance_metrics.total_sim_time = seconds(sim_end_time - sim_start_time);
%     performance_metrics.avg_iteration_time = mean(performance_metrics.sim_time);
%     performance_metrics.max_memory_usage_MB = max(performance_metrics.memory_usage_MB);
%     valid_cpu_readings = performance_metrics.cpu_usage(performance_metrics.cpu_usage >= 0);
%     if isempty(valid_cpu_readings), performance_metrics.avg_cpu_usage = 0;
%     else, performance_metrics.avg_cpu_usage = mean(valid_cpu_readings); end
% 
%     % 【修改】为输出文件名添加序号
%     output_filename = sprintf('simulation_performance_run_%d.txt', run_index);
%     try
%         fid = fopen(output_filename, 'w');
%         if fid == -1, error('无法打开文件 %s 进行写入。请检查文件权限。', output_filename); end
%         fprintf(fid, '=== Simulation Performance Metrics (Run %d) ===\n\n', run_index);
%         fprintf(fid, 'Total simulation time: %.2f seconds\n', performance_metrics.total_sim_time);
%         fprintf(fid, 'Average iteration time: %.4f ms\n', performance_metrics.avg_iteration_time * 1000);
%         fprintf(fid, 'Peak MATLAB process memory: %.2f MB\n', performance_metrics.max_memory_usage_MB);
%         fprintf(fid, 'Average system-wide CPU usage: %.2f%%\n\n', performance_metrics.avg_cpu_usage);
%         fprintf(fid, '=== Detailed Iteration Metrics ===\n');
%         fprintf(fid, '%-10s %-12s %-12s %-12s\n', 'Iteration', 'Time(ms)', 'Memory(MB)', 'CPU_Usage(%)');
%         fprintf(fid, '-----------------------------------------------------\n');
%         for i = 1:config.max_sim
%             fprintf(fid, '%-10d %-12.2f %-12.2f %-12.2f\n', i, ...
%                 performance_metrics.sim_time(i) * 1000, ...
%                 performance_metrics.memory_usage_MB(i), ...
%                 performance_metrics.cpu_usage(i));
%         end
%         fclose(fid);
%         fprintf('性能指标已成功保存到 %s\n', output_filename);
%     catch e
%         if fid ~= -1, fclose(fid); end
%         warning('写入性能文件时出错: %s', e.message);
%     end
% 
%     % 【新增改动】将本次运行的完整结果存入cell数组
%     all_run_results{run_index} = performance_metrics;
% 
%     %% 结果绘图
%     disp('正在生成结果图像...');
%     set(0,'defaultfigurecolor','w');
% 
%     % 【修改】为图像窗口添加序号，防止覆盖
%     figure('Name', sprintf('Co-visible Area of RS - Run %d', run_index));
%     plot(co_vis_area_stack, 'Linewidth', 2);
%     xlabel('Time(s)'); ylabel('Co-visible area'); title(sprintf('Area of co-visible area of RS (Run %d)', run_index)); grid on;
% 
%     figure('Name', sprintf('Payload Loss - Run %d', run_index));
%     hold on; 
%     for i=1:config.num_payload, plot(payload_loss_stack(i, :), 'Linewidth', 1.5); end
%     plot(ones([1,config.max_sim])*min(max(theta_max),0.5*pi - config.camera_view), 'Linewidth', 1.5, 'LineStyle', '--');
%     hold off;
%     xlabel('Time(s)'); ylabel('Payload loss(rad)'); title(sprintf('Angle between payload and camera direction of RS (Run %d)', run_index));
%     legendCell = cellstr(num2str((1:config.num_payload)', 'payload %d')); legendCell{end+1} = 'Safety Bound'; legend(legendCell, 'Location', 'best'); grid on;
% 
%     figure('Name', sprintf('Fuel Loss - Run %d', run_index));
%     hold on;
%     for i=2:2+config.num_rec_sat, plot(eny_consum_stack(i, :), 'Linewidth', 1.5); end
%     hold off;
%     xlabel('Time(s)'); ylabel('Fuel loss'); title(sprintf('Fuel loss of target satellite and reconnassiance satellite (Run %d)', run_index));
%     legend(["orb_TS", "orb_RS1", "orb_RS2", "orb_RS3"],'Interpreter','none', 'Location', 'best'); grid on;
% 
%     figure('Name', sprintf('Mission Loss - Run %d', run_index));
%     plot(mission_loss_stack, 'Linewidth', 1.5);
%     xlabel('Time(s)'); ylabel('mission loss(rad)'); title(sprintf('Angle between camera direction of TS and mission direction (Run %d)', run_index)); grid on;
% 
%     fprintf('第 %d / %d 次仿真运行完成。\n', run_index, num_runs);
% 
% % 【新增改动】结束N次运行的大循环
% end
% 
% disp('所有仿真运行全部完成！');
% 
% 
% %% ========================================================================
% %                     辅助模块：嵌套与本地函数 (无任何改动)
% % =========================================================================
% 
%     function update_cpu_usage_robust(~, ~)
%         usage = 0;
%         try
%             osBean = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
%             usage = osBean.getSystemCpuLoad() * 100;
%             if usage < 0, usage = 0; end
%         catch E
%             warning('MATLAB:CPUmonitor:JavaError', '通过Java获取CPU使用率失败 (%s)。将尝试使用系统命令。', E.message);
%             usage = fallback_cpu_usage();
%         end
%         if isempty(last_cpu_usages), last_cpu_usages = usage; end
%         current_cpu_usage = 0.7 * usage + 0.3 * last_cpu_usages;
%         last_cpu_usages = current_cpu_usage;
%     end
% 
% end 
% 
% 
% function stop_and_delete_timer(t)
%     if isvalid(t) && strcmp(t.Running, 'on'), stop(t); end
%     if isvalid(t), delete(t); disp('CPU monitor timer stopped and deleted.'); end
% end
% 
% function mem_bytes = get_process_resident_memory()
%     mem_bytes = 0;
%     try
%         if ispc
%             s = memory; mem_bytes = s.MemUsedMATLAB;
%         elseif isunix || ismac
%             pid = feature('getpid'); [~, result] = system(['ps -p ' num2str(pid) ' -o rss=']);
%             mem_kb = str2double(result); if ~isnan(mem_kb), mem_bytes = mem_kb * 1024; end
%         end
%     catch
%         mem_bytes = 0;
%     end
% end
% 
% function usage = fallback_cpu_usage()
%     usage = 0;
%     try
%         if ispc, [~,r] = system('wmic cpu get loadpercentage /value'); v = regexp(r,'LoadPercentage=(\d+)','tokens'); if ~isempty(v)&&~isempty(v{1}), usage=str2double(v{1}{1}); end
%         elseif isunix||ismac, [~,r] = system('top -bn1|grep "Cpu(s)"'); s = str2double(regexp(r,'[\d\.]+','match')); if ~isempty(s)&&length(s)>=4, usage=100-s(4); end
%         end
%     catch, usage = 0; end
%     if isempty(usage)||~isscalar(usage)||isnan(usage), usage = 0; end
% end