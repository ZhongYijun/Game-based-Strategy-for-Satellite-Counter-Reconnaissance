% % =========================================================================
% %           FINAL MAIN SCRIPT FOR ALGORITHM COMPARISON (REVISED)
% % =========================================================================
% % This script compares two re-implemented controller philosophies in a
% % ONE-ON-ONE game scenario, with all identified bugs fixed.
% % 1. Dong et al. Style (GA)
% % 2. Wang et al. Style (PSO)
% % =========================================================================
% 
% clear; clc; close all;
% fprintf('=== ALGORITHM COMPARISON (FIXED): Dong (GA) vs. Wang (PSO) ===\n');
% 
% %% 1. CONFIGURATION
% config = PP_Config();
% if config.num_rec_sat ~= 1
%     warning('Overriding config: Setting num_rec_sat to 1 for this 1-on-1 comparison.');
%     config.num_rec_sat = 1;
% end
% % --- Parameters to control EA solver computational cost PER DECISION STEP ---
% control_params.ga_pop_size = 100;
% control_params.ga_max_gen = 50;
% control_params.pso_swarm_size = 100;
% control_params.pso_max_iter = 50;
% 
% %% 2. INITIALIZATION FOR TWO PARALLEL SIMULATIONS
% initial_state = [config.tgt_att_0; config.rel_orb_0(:,1)]; 
% 
% % --- Storages for Dong et al. ---
% state_hist_dong = zeros(size(initial_state, 1), config.max_sim + 1);
% state_hist_dong(:, 1) = initial_state;
% co_vis_area_hist_dong = zeros(1, config.max_sim);
% payload_loss_hist_dong = pi*ones(config.num_payload, config.max_sim);
% fuel_hist_dong = zeros(2, config.max_sim); % 1 for Target Orb, 1 for Recon Orb
% mission_loss_hist_dong = zeros(1, config.max_sim);
% time_hist_dong = zeros(1, config.max_sim);
% 
% % --- Storages for Wang et al. ---
% state_hist_wang = zeros(size(initial_state, 1), config.max_sim + 1);
% state_hist_wang(:, 1) = initial_state;
% co_vis_area_hist_wang = zeros(1, config.max_sim);
% payload_loss_hist_wang = pi*ones(config.num_payload, config.max_sim);
% fuel_hist_wang = zeros(2, config.max_sim);
% mission_loss_hist_wang = zeros(1, config.max_sim);
% time_hist_wang = zeros(1, config.max_sim);
% 
% % Shared variables
% unit_ball_sample = unit_ball_sample_gen(config);
% non_zero_count_dong = zeros(config.max_num_ball_sample, 1);
% non_zero_count_wang = zeros(config.max_num_ball_sample, 1);
% 
% %% 3. MAIN SIMULATION LOOP
% fprintf('\n--- Starting Parallel Simulations for %d steps ---\n', config.max_sim);
% for i = 1:config.max_sim
%     if mod(i, 100) == 0, fprintf('--- Step %d / %d ---\n', i, config.max_sim); end
% 
%     theta = config.orb_ele(6) + i*config.sample_time_span*config.orb_w; 
%     J2L_mat = calc_j2l_mat(config, theta);
% 
%     % --- RUN CONTROLLER 1: DONG ---
%     tic;
%     control_dong = controller_Dong(config, state_hist_dong(:, i), J2L_mat, unit_ball_sample, non_zero_count_dong, control_params);
%     time_hist_dong(i) = toc;
% 
%     % --- RUN CONTROLLER 2: WANG ---
%     tic;
%     control_wang = controller_Wang(config, state_hist_wang(:, i), J2L_mat, unit_ball_sample, non_zero_count_wang, control_params);
%     time_hist_wang(i) = toc;
% 
%     % --- SIMULATE & LOG ---
%     state_hist_dong(:, i+1) = Simulator(config, state_hist_dong(:, i), control_dong);
%     state_hist_wang(:, i+1) = Simulator(config, state_hist_wang(:, i), control_wang);
% 
%     [co_vis_area_hist_dong(i), payload_loss_hist_dong(:,i), fuel_hist_dong(:,i), mission_loss_hist_dong(i), non_zero_count_dong] = ...
%         log_metrics(config, state_hist_dong(:, i), control_dong, J2L_mat, unit_ball_sample, non_zero_count_dong);
% 
%     [co_vis_area_hist_wang(i), payload_loss_hist_wang(:,i), fuel_hist_wang(:,i), mission_loss_hist_wang(i), non_zero_count_wang] = ...
%         log_metrics(config, state_hist_wang(:, i), control_wang, J2L_mat, unit_ball_sample, non_zero_count_wang);
% end
% fprintf('--- All simulations finished. ---\n');
% 
% %% 4. PLOT COMPARISON RESULTS
% time_axis = 1:config.max_sim;
% set(0,'defaultfigurecolor','w');
% 
% figure('Name', 'Comparison: Co-visible Area'); hold on; plot(time_axis, co_vis_area_hist_dong, 'r--', 'LineWidth', 1.5); plot(time_axis, co_vis_area_hist_wang, 'g-.', 'LineWidth', 1.5); hold off; title('Co-visible Area'); xlabel('Time (s)'); ylabel('Area (steradians)'); legend('Dong et al. (GA)', 'Wang et al. (PSO)'); grid on;
% figure('Name', 'Comparison: Payload Loss'); hold on; plot(time_axis, payload_loss_hist_dong(1, :), 'r--', 'LineWidth', 1.5); plot(time_axis, payload_loss_hist_wang(1, :), 'g-.', 'LineWidth', 1.5); theta_max_final = safe_acos(config.radius_tgt_sat/norm(state_hist_dong(config.att_dim+1:config.att_dim+3, end))); plot(time_axis, ones(1,config.max_sim)*min(theta_max_final, 0.5*pi - config.camera_view), 'k:', 'LineWidth', 1); hold off; title('Payload Loss Comparison (Payload 1)'); xlabel('Time (s)'); ylabel('Loss (rad)'); legend('Dong (GA)', 'Wang (PSO)', 'Safety Bound'); grid on;
% figure('Name', 'Comparison: Cumulative Fuel Consumption'); hold on; plot(time_axis, cumsum(fuel_hist_dong(1,:)), 'r--','DisplayName','Dong (GA) - Target'); plot(time_axis, cumsum(fuel_hist_dong(2,:)), 'r-','DisplayName','Dong (GA) - Recon'); plot(time_axis, cumsum(fuel_hist_wang(1,:)), 'g--','DisplayName','Wang (PSO) - Target'); plot(time_axis, cumsum(fuel_hist_wang(2,:)), 'g-','DisplayName','Wang (PSO) - Recon'); hold off; title('Cumulative Orbital Fuel Consumption (\DeltaV)'); xlabel('Time (s)'); ylabel('Cumulative \DeltaV (m/s)'); legend; grid on;
% figure('Name', 'Comparison: Mission Loss'); hold on; plot(time_axis, mission_loss_hist_dong, 'r--', 'LineWidth', 1.5); plot(time_axis, mission_loss_hist_wang, 'g-.', 'LineWidth', 1.5); hold off; title('Mission Loss of Target Satellite'); xlabel('Time (s)'); ylabel('Loss (rad)'); legend('Dong et al. (GA)', 'Wang et al. (PSO)'); grid on;
% figure('Name', 'Comparison: Computation Time per Step'); hold on; plot(time_axis, time_hist_dong * 1000, 'r--', 'LineWidth', 1.5); plot(time_axis, time_hist_wang * 1000, 'g-.', 'LineWidth', 1.5); hold off; avg_time_dong_ms = mean(time_hist_dong) * 1000; avg_time_wang_ms = mean(time_hist_wang) * 1000; title(sprintf('Computation Time per Step (Avg: Dong %.1f ms, Wang %.1f ms)', avg_time_dong_ms, avg_time_wang_ms)); xlabel('Time (s)'); ylabel('Time per Step (ms)'); legend('Dong et al. (GA)', 'Wang et al. (PSO)'); grid on; set(gca, 'YScale', 'log');
% 
% %% Helper Functions
% function J2L = calc_j2l_mat(config, theta), i=config.orb_ele(3); O=config.orb_ele(4); ct=cos(theta);st=sin(theta);ci=cos(i);si=sin(i);cO=cos(O);sO=sin(O); J2L=[ct*cO-st*sO*ci,ct*sO+st*cO*ci,st*si;-st*cO-ct*sO*ci,-st*sO+ct*cO*ci,ct*si;si*sO,-si*cO,ci]; end
% function [co,pl,fuel,ml,nzc_out] = log_metrics(c,x,u,J2L,usb,nzc_in), [tm,cr]=att_trans_mat(c,x); nzc_out=non_zero_count_updater(c,x,cr,J2L,usb,nzc_in); co=4*pi/c.max_num_ball_sample*sum(nzc_out); pl=pi*ones(c.num_payload,1); tm_log=safe_acos(c.radius_tgt_sat/norm(x(c.att_dim+1:c.att_dim+3))); for j=1:c.num_payload,if cr'*J2L*c.a_sun<=cos(c.ang_sun),pl(j)=min(pl(j),safe_acos(c.q(:,j)'*tm*cr));end,end; fuel=[norm(u(c.ctl_dim+1:2*c.ctl_dim));norm(u(2*c.ctl_dim+1:end))]; ml=safe_acos(c.q(:,1)'*tm*c.a_mission); end
% function nzc_up=non_zero_count_updater(c,x,cr,J2L,usb,nzc_in), nzc_up=nzc_in;tm_upd=safe_acos(c.radius_tgt_sat/norm(x(c.att_dim+1:c.att_dim+3))); if norm(x(c.att_dim+1:c.att_dim+3))<=c.eff_obsv_distance, is_gv=(usb*cr>=cos(min(tm_upd,0.5*pi-c.camera_view))); is_io=(usb*J2L*c.a_sun<=cos(c.ang_sun)); nzc_up=nzc_up|(is_gv&is_io); end, end
% function angle=safe_acos(x), angle=acos(max(-1.0,min(1.0,x))); end

% =========================================================================
%           FINAL MAIN SCRIPT FOR ALGORITHM COMPARISON (3-WAY FIXED)
% =========================================================================
% This script systematically compares THREE different controller philosophies
% in a unified 1-on-1 game scenario, with all identified bugs and logical flaws corrected.
% 1. Our Original MPC    : The fmincon-based controller from the paper.
% 2. Dong et al. Style   : A time-driven approach using a Genetic Algorithm.
% 3. Wang et al. Style   : A utility-driven approach using Particle Swarm Opt.
% =========================================================================

clear; clc; close all;
fprintf('=== ALGORITHM COMPARISON (3-WAY FINAL): Ours (MPC) vs. Dong (GA) vs. Wang (PSO) ===\n');

%% 1. CONFIGURATION
% =========================================================================
config = PP_Config();
if config.num_rec_sat ~= 1
    warning('Overriding config: Setting num_rec_sat to 1 for this 1-on-1 comparison.');
    config.num_rec_sat = 1;
end
% --- Parameters to control EA solver computational cost ---
control_params.ga_pop_size = 100;
control_params.ga_max_gen = 50;
control_params.pso_swarm_size = 100;
control_params.pso_max_iter = 50;

%% 2. INITIALIZATION FOR THREE PARALLEL SIMULATIONS
% =========================================================================
initial_state = [config.tgt_att_0; config.rel_orb_0(:,1)]; 

% --- Storages for "Our MPC Controller" ---
state_hist_ours = zeros(size(initial_state, 1), config.max_sim + 1); state_hist_ours(:, 1) = initial_state;
co_vis_area_hist_ours = zeros(1, config.max_sim); payload_loss_hist_ours = pi*ones(config.num_payload, config.max_sim);
fuel_hist_ours = zeros(2, config.max_sim); mission_loss_hist_ours = zeros(1, config.max_sim); time_hist_ours = zeros(1, config.max_sim);
cur_ctl_ours = zeros([(2+config.num_rec_sat)*config.ctl_dim, config.num_sample+1]);
state_sample_stack_ours = zeros(size(initial_state, 1), config.num_sample);

% --- Storages for "Dong et al. Style Controller" ---
state_hist_dong = zeros(size(initial_state, 1), config.max_sim + 1); state_hist_dong(:, 1) = initial_state;
co_vis_area_hist_dong = zeros(1, config.max_sim); payload_loss_hist_dong = pi*ones(config.num_payload, config.max_sim);
fuel_hist_dong = zeros(2, config.max_sim); mission_loss_hist_dong = zeros(1, config.max_sim); time_hist_dong = zeros(1, config.max_sim);

% --- Storages for "Wang et al. Style Controller" ---
state_hist_wang = zeros(size(initial_state, 1), config.max_sim + 1); state_hist_wang(:, 1) = initial_state;
co_vis_area_hist_wang = zeros(1, config.max_sim); payload_loss_hist_wang = pi*ones(config.num_payload, config.max_sim);
fuel_hist_wang = zeros(2, config.max_sim); mission_loss_hist_wang = zeros(1, config.max_sim); time_hist_wang = zeros(1, config.max_sim);

% Shared variables
unit_ball_sample = unit_ball_sample_gen(config);
non_zero_count_ours = zeros(config.max_num_ball_sample, 1);
non_zero_count_dong = zeros(config.max_num_ball_sample, 1);
non_zero_count_wang = zeros(config.max_num_ball_sample, 1);

%% 3. MAIN SIMULATION LOOP
% =========================================================================
fprintf('\n--- Starting Parallel Simulations for %d steps ---\n', config.max_sim);
for i = 1:config.max_sim
    if mod(i, 100) == 0, fprintf('--- Step %d / %d ---\n', i, config.max_sim); end

    theta = config.orb_ele(6) + i*config.sample_time_span*config.orb_w; 
    J2L_mat = calc_j2l_mat(config, theta);
    
    % --- STATE UPDATE LOGIC: Update history for EACH simulation branch BEFORE decisions ---
    [~, cur_r_ours] = att_trans_mat(config, state_hist_ours(:, i));
    non_zero_count_ours = non_zero_count_updater(config, state_hist_ours(:, i), cur_r_ours, J2L_mat, unit_ball_sample, non_zero_count_ours);
    
    [~, cur_r_dong] = att_trans_mat(config, state_hist_dong(:, i));
    non_zero_count_dong = non_zero_count_updater(config, state_hist_dong(:, i), cur_r_dong, J2L_mat, unit_ball_sample, non_zero_count_dong);
    
    [~, cur_r_wang] = att_trans_mat(config, state_hist_wang(:, i));
    non_zero_count_wang = non_zero_count_updater(config, state_hist_wang(:, i), cur_r_wang, J2L_mat, unit_ball_sample, non_zero_count_wang);
    
    % --- RUN CONTROLLERS ---
    tic; cur_ctl_ours = mpc_controller(config, state_hist_ours(:, i), state_sample_stack_ours, cur_ctl_ours, J2L_mat, unit_ball_sample, non_zero_count_ours); time_hist_ours(i) = toc;
    control_to_apply_ours = cur_ctl_ours(:,1);

    tic; control_dong = controller_Dong(config, state_hist_dong(:, i), J2L_mat, unit_ball_sample, non_zero_count_dong, control_params); time_hist_dong(i) = toc;
    
    tic; control_wang = controller_Wang(config, state_hist_wang(:, i), J2L_mat, unit_ball_sample, non_zero_count_wang, control_params); time_hist_wang(i) = toc;
    
    % --- SIMULATE NEXT STATE ---
    state_hist_ours(:, i+1) = Simulator(config, state_hist_ours(:, i), control_to_apply_ours);
    state_hist_dong(:, i+1) = Simulator(config, state_hist_dong(:, i), control_dong);
    state_hist_wang(:, i+1) = Simulator(config, state_hist_wang(:, i), control_wang);

    % --- LOG METRICS ---
    [co_vis_area_hist_ours(i), payload_loss_hist_ours(:,i), fuel_hist_ours(:,i), mission_loss_hist_ours(i)] = log_metrics(config, state_hist_ours(:, i), control_to_apply_ours, J2L_mat, non_zero_count_ours);
    [co_vis_area_hist_dong(i), payload_loss_hist_dong(:,i), fuel_hist_dong(:,i), mission_loss_hist_dong(i)] = log_metrics(config, state_hist_dong(:, i), control_dong, J2L_mat, non_zero_count_dong);
    [co_vis_area_hist_wang(i), payload_loss_hist_wang(:,i), fuel_hist_wang(:,i), mission_loss_hist_wang(i)] = log_metrics(config, state_hist_wang(:, i), control_wang, J2L_mat, non_zero_count_wang);
end
fprintf('--- All simulations finished. ---\n');

%% 4. PLOT COMPARISON RESULTS
% =========================================================================
time_axis = 1:config.max_sim; set(0,'defaultfigurecolor','w');
figure('Name','Comparison: Co-visible Area'); hold on; plot(time_axis,co_vis_area_hist_ours,'b-','LineWidth',2); plot(time_axis,co_vis_area_hist_dong,'r--','LineWidth',1.5); plot(time_axis,co_vis_area_hist_wang,'g-.','LineWidth',1.5); hold off; title('Co-visible Area'); xlabel('Time (s)'); ylabel('Area (steradians)'); legend('Our MPC','Dong (GA)','Wang (PSO)'); grid on;
figure('Name','Comparison: Payload Loss'); hold on; plot(time_axis,payload_loss_hist_ours(1,:),'b-','LineWidth',2); plot(time_axis,payload_loss_hist_dong(1,:),'r--','LineWidth',1.5); plot(time_axis,payload_loss_hist_wang(1,:),'g-.','LineWidth',1.5); theta_max_final=safe_acos(config.radius_tgt_sat/norm(state_hist_ours(config.att_dim+1:config.att_dim+3,end))); plot(time_axis,ones(1,config.max_sim)*min(theta_max_final,0.5*pi-config.camera_view),'k:','LineWidth',1); hold off; title('Payload Loss (Payload 1)'); xlabel('Time (s)'); ylabel('Loss (rad)'); legend('Our MPC','Dong (GA)','Wang (PSO)','Safety Bound'); grid on;
figure('Name','Comparison: Cumulative Fuel Consumption'); hold on; plot(time_axis,cumsum(fuel_hist_ours(1,:)),'b--','DisplayName','Our MPC - Target'); plot(time_axis,cumsum(fuel_hist_ours(2,:)),'b-','DisplayName','Our MPC - Recon'); plot(time_axis,cumsum(fuel_hist_dong(1,:)),'r--','DisplayName','Dong (GA) - Target'); plot(time_axis,cumsum(fuel_hist_dong(2,:)),'r-','DisplayName','Dong (GA) - Recon'); plot(time_axis,cumsum(fuel_hist_wang(1,:)),'g--','DisplayName','Wang (PSO) - Target'); plot(time_axis,cumsum(fuel_hist_wang(2,:)),'g-','DisplayName','Wang (PSO) - Recon'); hold off; title('Cumulative Orbital Fuel Consumption'); xlabel('Time (s)'); ylabel('Cumulative \DeltaV (m/s)'); legend('Location','northwest'); grid on;
figure('Name','Comparison: Mission Loss'); hold on; plot(time_axis,mission_loss_hist_ours,'b-','LineWidth',2); plot(time_axis,mission_loss_hist_dong,'r--','LineWidth',1.5); plot(time_axis,mission_loss_hist_wang,'g-.','LineWidth',1.5); hold off; title('Mission Loss'); xlabel('Time (s)'); ylabel('Loss (rad)'); legend('Our MPC','Dong (GA)','Wang (PSO)'); grid on;
figure('Name','Comparison: Computation Time'); hold on; plot(time_axis,time_hist_ours*1000,'b-','LineWidth',2); plot(time_axis,time_hist_dong*1000,'r--','LineWidth',1.5); plot(time_axis,time_hist_wang*1000,'g-.','LineWidth',1.5); hold off; avg_ours=mean(time_hist_ours)*1000; avg_dong=mean(time_hist_dong)*1000; avg_wang=mean(time_hist_wang)*1000; title(sprintf('Computation Time (Avg: Ours %.1f, Dong %.1f, Wang %.1f ms)',avg_ours,avg_dong,avg_wang)); xlabel('Time (s)'); ylabel('Time per Step (ms)'); legend('Our MPC','Dong (GA)','Wang (PSO)'); grid on; set(gca,'YScale','log');

%% Helper Functions
function J2L=calc_j2l_mat(c,t),i=c.orb_ele(3);O=c.orb_ele(4);ct=cos(t);st=sin(t);ci=cos(i);si=sin(i);cO=cos(O);sO=sin(O);J2L=[ct*cO-st*sO*ci,ct*sO+st*cO*ci,st*si;-st*cO-ct*sO*ci,-st*sO+ct*cO*ci,ct*si;si*sO,-si*cO,ci];end
function [co,pl,fuel,ml]=log_metrics(c,x,u,J2L,nzc),[tm,cr]=att_trans_mat(c,x);co=4*pi/c.max_num_ball_sample*sum(nzc);pl=pi*ones(c.num_payload,1);rel_norm=norm(x(c.att_dim+1:c.att_dim+3));if rel_norm<1e-6,rel_norm=1e-6;end;tm_log=safe_acos(c.radius_tgt_sat/rel_norm);for j=1:c.num_payload,if cr'*J2L*c.a_sun<=cos(c.ang_sun),pl(j)=min(pl(j),safe_acos(c.q(:,j)'*tm*cr));end,end;fuel=[norm(u(c.ctl_dim+1:2*c.ctl_dim));norm(u(2*c.ctl_dim+1:end))];ml=safe_acos(c.q(:,1)'*tm*c.a_mission);end
function nzc_up=non_zero_count_updater(c,x,cr,J2L,usb,nzc_in),nzc_up=nzc_in;rel_pos=x(c.att_dim+1:c.att_dim+3);rel_norm=norm(rel_pos);if rel_norm<1e-6,return;end;tm_upd=safe_acos(c.radius_tgt_sat/rel_norm);if rel_norm<=c.eff_obsv_distance,is_gv=(usb*cr>=cos(min(tm_upd,0.5*pi-c.camera_view)));is_io=(usb*J2L*c.a_sun<=cos(c.ang_sun));nzc_up=nzc_up|(is_gv&is_io);end,end
function angle=safe_acos(x),angle=acos(max(-1.0,min(1.0,x)));end