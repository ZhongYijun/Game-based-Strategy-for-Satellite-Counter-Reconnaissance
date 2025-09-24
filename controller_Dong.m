% Filename: controller_Dong.m (CORRECTED & FINAL VERSION)

function [control_vector] = controller_Dong(config, cur_state, J2L_mat, unit_ball_sample, non_zero_count, control_params)
% CONTROLLER: Implements the "Time-Driven" philosophy from Dong et al.
% USES: Genetic Algorithm (GA) to solve an ORBITAL-ONLY game in a 1-on-1 scenario.
% Attitude control is always ZERO. Norm constraints and history are correctly handled.

if nargin < 6, control_params.ga_pop_size=50; control_params.ga_max_gen=20; else
    if ~isfield(control_params,'ga_pop_size'),control_params.ga_pop_size=50;end
    if ~isfield(control_params,'ga_max_gen'),control_params.ga_max_gen=20;end
end
horizon=config.T_sampling;control_dim_orb=2*config.ctl_dim;total_optim_vars=horizon*control_dim_orb;
lb=-ones(total_optim_vars,1);ub=ones(total_optim_vars,1);
opts=optimoptions('ga','PopulationSize',control_params.ga_pop_size,'MaxGenerations',control_params.ga_max_gen,'Display','off','UseParallel',true);
fitness_handle = @(x) fitness_function_Dong(x, config, cur_state, non_zero_count, unit_ball_sample);
[solution_vector,~]=ga(fitness_handle,total_optim_vars,[],[],[],[],lb,ub,[],opts);
optimal_orb_sequence=reshape(solution_vector,control_dim_orb,horizon);
u_orb_all_norm=optimal_orb_sequence(:,1);
u_orb_tgt_raw=u_orb_all_norm(1:3)*config.tgt_orb_u_max;u_orb_rec_raw=u_orb_all_norm(4:6)*config.rec_orb_u_max(1);
if norm(u_orb_tgt_raw)>config.tgt_orb_u_max,u_orb_tgt=u_orb_tgt_raw/norm(u_orb_tgt_raw)*config.tgt_orb_u_max;else,u_orb_tgt=u_orb_tgt_raw;end
if norm(u_orb_rec_raw)>config.rec_orb_u_max(1),u_orb_rec=u_orb_rec_raw/norm(u_orb_rec_raw)*config.rec_orb_u_max(1);else,u_orb_rec=u_orb_rec_raw;end
control_vector=[zeros(config.ctl_dim,1);u_orb_tgt;u_orb_rec];
end

function total_cost = fitness_function_Dong(orb_control_vector, config, start_state, start_non_zero_count, unit_ball_sample)
    control_dim_orb=2*config.ctl_dim;orb_control_sequence_norm=reshape(orb_control_vector,control_dim_orb,config.T_sampling);
    total_fuel=0;cur_state=start_state;nzc=start_non_zero_count;
    for t=1:config.T_sampling
        u_att_zero=zeros(config.ctl_dim,1);cur_orb_ctl_norm=orb_control_sequence_norm(:,t);
        u_orb_tgt_raw=cur_orb_ctl_norm(1:3)*config.tgt_orb_u_max;u_orb_rec_raw=cur_orb_ctl_norm(4:6)*config.rec_orb_u_max(1);
        if norm(u_orb_tgt_raw)>config.tgt_orb_u_max,u_orb_tgt=u_orb_tgt_raw/norm(u_orb_tgt_raw)*config.tgt_orb_u_max;else,u_orb_tgt=u_orb_tgt_raw;end
        if norm(u_orb_rec_raw)>config.rec_orb_u_max(1),u_orb_rec=u_orb_rec_raw/norm(u_orb_rec_raw)*config.rec_orb_u_max(1);else,u_orb_rec=u_orb_rec_raw;end
        cur_ctl_t=[u_att_zero;u_orb_tgt;u_orb_rec];
        
        total_fuel=total_fuel+norm(u_orb_tgt)+norm(u_orb_rec);
        cur_state=Simulator(config,cur_state,cur_ctl_t);
        [~,cr]=att_trans_mat(config,cur_state);J2L_mat=calc_j2l_mat_local(config,t);
        nzc=non_zero_count_updater_local(config,cur_state,cr,J2L_mat,unit_ball_sample,nzc);
    end
    final_co_vis_area=4*pi/config.max_num_ball_sample*sum(nzc);
    total_cost=1e4*final_co_vis_area+1.0*total_fuel;
end
function nzc_up=non_zero_count_updater_local(c,x,cr,J2L,usb,nzc_in),nzc_up=nzc_in;rel_pos=x(c.att_dim+1:c.att_dim+3);rel_norm=norm(rel_pos);if rel_norm<1e-6,return;end;tm_upd=safe_acos(c.radius_tgt_sat/rel_norm);if rel_norm<=c.eff_obsv_distance,is_gv=(usb*cr>=cos(min(tm_upd,0.5*pi-c.camera_view)));is_io=(usb*J2L*c.a_sun<=cos(c.ang_sun));nzc_up=nzc_up|(is_gv&is_io);end,end
function J2L=calc_j2l_mat_local(c,t_off),theta=c.orb_ele(6)+t_off*c.sample_time_span*c.orb_w;i=c.orb_ele(3);O=c.orb_ele(4);ct=cos(theta);st=sin(theta);ci=cos(i);si=sin(i);cO=cos(O);sO=sin(O);J2L=[ct*cO-st*sO*ci,ct*sO+st*cO*ci,st*si;-st*cO-ct*sO*ci,-st*sO+ct*cO*ci,ct*si;si*sO,-si*cO,ci];end
function a=safe_acos(x),a=acos(max(-1,min(1,x)));end