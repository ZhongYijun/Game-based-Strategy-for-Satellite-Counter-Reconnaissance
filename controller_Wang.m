% Filename: controller_Wang.m (CORRECTED & FINAL VERSION)

function [control_vector] = controller_Wang(config, cur_state, J2L_mat, unit_ball_sample, non_zero_count, control_params)
% CONTROLLER: Implements the "Utility-Driven" philosophy from Wang et al.
% USES: Particle Swarm Opt. (PSO) to solve an ORBITAL-ONLY game in a 1-on-1 scenario.
% Attitude control is always ZERO. Norm constraints and history are correctly handled.

if nargin < 6, control_params.pso_swarm_size=50; control_params.pso_max_iter=30; else
    if ~isfield(control_params,'pso_swarm_size'),control_params.pso_swarm_size=50;end
    if ~isfield(control_params,'pso_max_iter'),control_params.pso_max_iter=30;end
end
horizon=config.T_sampling;control_dim_orb=2*config.ctl_dim;total_optim_vars=horizon*control_dim_orb;
lb=-ones(total_optim_vars,1);ub=ones(total_optim_vars,1);
opts=optimoptions('particleswarm','SwarmSize',control_params.pso_swarm_size,'MaxIterations',control_params.pso_max_iter,'Display','off','UseParallel',true);
fitness_handle = @(x) fitness_function_Wang(x, config, cur_state);
[solution_vector,~]=particleswarm(fitness_handle,total_optim_vars,lb,ub,opts);
optimal_orb_sequence=reshape(solution_vector,control_dim_orb,horizon);
u_orb_all_norm=optimal_orb_sequence(:,1);
u_orb_tgt_raw=u_orb_all_norm(1:3)*config.tgt_orb_u_max;u_orb_rec_raw=u_orb_all_norm(4:6)*config.rec_orb_u_max(1);
if norm(u_orb_tgt_raw)>config.tgt_orb_u_max,u_orb_tgt=u_orb_tgt_raw/norm(u_orb_tgt_raw)*config.tgt_orb_u_max;else,u_orb_tgt=u_orb_tgt_raw;end
if norm(u_orb_rec_raw)>config.rec_orb_u_max(1),u_orb_rec=u_orb_rec_raw/norm(u_orb_rec_raw)*config.rec_orb_u_max(1);else,u_orb_rec=u_orb_rec_raw;end
control_vector=[zeros(config.ctl_dim,1);u_orb_tgt;u_orb_rec];
end

function total_cost = fitness_function_Wang(orb_control_vector, config, start_state)
    % Fitness function for the Utility-Driven (Wang et al.) philosophy.
    control_dim_orb=2*config.ctl_dim;orb_control_sequence_norm=reshape(orb_control_vector,control_dim_orb,config.T_sampling);
    total_com_utility=0;total_fuel=0;cur_state=start_state;
    for t=1:config.T_sampling
        u_att_zero=zeros(config.ctl_dim,1);cur_orb_ctl_norm=orb_control_sequence_norm(:,t);
        u_orb_tgt_raw=cur_orb_ctl_norm(1:3)*config.tgt_orb_u_max;u_orb_rec_raw=cur_orb_ctl_norm(4:6)*config.rec_orb_u_max(1);
        if norm(u_orb_tgt_raw)>config.tgt_orb_u_max,u_orb_tgt=u_orb_tgt_raw/norm(u_orb_tgt_raw)*config.tgt_orb_u_max;else,u_orb_tgt=u_orb_tgt_raw;end
        if norm(u_orb_rec_raw)>config.rec_orb_u_max(1),u_orb_rec=u_orb_rec_raw/norm(u_orb_rec_raw)*config.rec_orb_u_max(1);else,u_orb_rec=u_orb_rec_raw;end
        cur_ctl_t=[u_att_zero;u_orb_tgt;u_orb_rec];
        
        rel_pos_CoM=cur_state(config.att_dim+1:config.att_dim+3);dist_CoM=norm(rel_pos_CoM);
        if dist_CoM<1e-6,total_cost=1e10;return;end
        
        f_dist=1/(1+abs(dist_CoM-config.eff_obsv_distance));
        J2L_mat=calc_j2l_mat_local(config,t);sun_vec_lvlh=J2L_mat*config.a_sun;
        cos_theta_illum=dot(sun_vec_lvlh,-rel_pos_CoM/dist_CoM);
        f_illum=(1+cos_theta_illum)/2;
        utility_this_step=f_dist+f_illum;
        total_com_utility=total_com_utility+utility_this_step;
        total_fuel=total_fuel+norm(u_orb_tgt)+norm(u_orb_rec);
        cur_state=Simulator(config,cur_state,cur_ctl_t);
    end
    total_cost=100*total_com_utility+10*total_fuel;
end
function J2L=calc_j2l_mat_local(c,t_off),theta=c.orb_ele(6)+t_off*c.sample_time_span*c.orb_w;i=c.orb_ele(3);O=c.orb_ele(4);ct=cos(theta);st=sin(theta);ci=cos(i);si=sin(i);cO=cos(O);sO=sin(O);J2L=[ct*cO-st*sO*ci,ct*sO+st*cO*ci,st*si;-st*cO-ct*sO*ci,-st*sO+ct*cO*ci,ct*si;si*sO,-si*cO,ci];end