% function [ini_admissible_ctl] = mpc_controller(config, cur_state,...
%     state_sample_stack, cur_ctl, J2L_mat, unit_ball_sample, non_zero_count)
% 
% 
% ini_admissible_ctl = cur_ctl;
% 
% theta_max = zeros([config.num_rec_sat, 1]);
% for j=1:config.num_rec_sat
%     theta_max(j,1) = acos(config.radius_tgt_sat/norm(cur_state(...
%            config.att_dim+(j-1)*config.orb_dim+(1:3), 1)));
% end
% 
% opts = optimoptions("fmincon", 'Algorithm', "sqp", 'MaxIterations', 600,'MaxFunctionEvaluations',600);
% 
% for i=1:config.num_sample + 1
% 
%     ini_ctl_seq = repmat(cur_ctl(:,i), 1, config.T_sampling);
%     if i==1
% 
%         ctl_seq_r = fmincon(@(ctl_seq) loss_fun(config, cur_state, ...
%             [ini_ctl_seq(1:2*config.ctl_dim, :);ctl_seq], J2L_mat, ...
%              "V_R",unit_ball_sample, non_zero_count), ini_ctl_seq( ...
%              2*config.ctl_dim+1:end, :),[],[],[],[],[],[], @(ctl_seq) ...
%              mpc_cst(config, ctl_seq, "V_R"), opts);
%         ctl_seq_r = zeros([config.num_rec_sat*config.ctl_dim,config.T_sampling]);
% 
%         ctl_seq_t = fmincon(@(ctl_seq) loss_fun(config, cur_state, ...
%             [ctl_seq; ctl_seq_r], J2L_mat, "V_T", unit_ball_sample, ...
%             non_zero_count), ini_ctl_seq(1:2*config.ctl_dim, :), [],[],[], ...
%             [],[],[], @(ctl_seq) mpc_cst(config, ctl_seq, "V_T"), opts);
% 
%         ini_admissible_ctl(:, 1) = [ctl_seq_t(:,1);ctl_seq_r(:,1)];
%     else
% 
%         ctl_seq_r = fmincon(@(ctl_seq) loss_fun(config, state_sample_stack(...
%             :, i-1), [ini_ctl_seq(1:2*config.ctl_dim, :);ctl_seq], ...
%             J2L_mat, "V_R", unit_ball_sample, non_zero_count), ini_ctl_seq( ...
%             2*config.ctl_dim+1: end, :), [],[],[],[],[],[], @(ctl_seq) ...
%             mpc_cst(config, ctl_seq, "V_R"), opts);
%         ctl_seq_r = zeros([config.num_rec_sat*config.ctl_dim,config.T_sampling]);
% 
%         ctl_seq_t = fmincon(@(ctl_seq) loss_fun(config, state_sample_stack(...
%             :, i-1), [ctl_seq; ctl_seq_r], J2L_mat, "V_T", ...
%             unit_ball_sample, non_zero_count), ini_ctl_seq(1: 2*config.ctl_dim, ...
%             :), [],[],[],[],[],[], @(ctl_seq) mpc_cst(config, ctl_seq, "V_T"), opts);
% 
%         ini_admissible_ctl(:, i) = [ctl_seq_t(:,1); ctl_seq_r(:,1)];
%     end
%     if norm(ini_admissible_ctl(1:config.ctl_dim, i)) > config.tgt_att_u_max
% 
%             ini_admissible_ctl(1:config.ctl_dim, i) = config.tgt_att_u_max*...
%             ini_admissible_ctl(1:config.ctl_dim, i)/norm(ini_admissible_ctl...
%             (1:config.ctl_dim, i));
%     end
% 
%     norm_u_max = config.tgt_orb_u_max;
%     for k = 1:config.num_rec_sat+1
%         if norm(ini_admissible_ctl(config.ctl_dim*k + (1:config.ctl_dim), i)) > ...
%                 config.tgt_orb_u_max
%             if k > 1
%                 norm_u_max = config.rec_orb_u_max(k-1);
%             end
%             ini_admissible_ctl(config.ctl_dim*k + (1:config.ctl_dim), i) = ...
%                 norm_u_max*ini_admissible_ctl(config.ctl_dim*k +...
%                 (1:config.ctl_dim), i)/norm(ini_admissible_ctl(config.ctl_dim*k...
%                 + (1:config.ctl_dim), i));
%         end
%     end
% 
% end
% 
% end

% Filename: mpc_controller.m (REVISED to support EA solvers)
% You can save this as a new file, e.g., mpc_controller_EA.m

function [ini_admissible_ctl] = mpc_controller(config, cur_state,...
    state_sample_stack, cur_ctl, J2L_mat, unit_ball_sample, non_zero_count, solver_type)

if nargin < 8
    solver_type = 'fmincon'; % Default to your original solver
end

ini_admissible_ctl = cur_ctl;

theta_max = zeros([config.num_rec_sat, 1]);
for j=1:config.num_rec_sat
    theta_max(j,1) = acos(config.radius_tgt_sat/norm(cur_state(...
           config.att_dim+(j-1)*config.orb_dim+(1:3), 1)));
end

% --- START OF REFACTORING ---

% Set options for different solvers
if strcmpi(solver_type, 'fmincon')
    opts = optimoptions("fmincon", 'Algorithm', "sqp", 'MaxIterations', 500, 'MaxFunctionEvaluations', 500);
elseif strcmpi(solver_type, 'ga')
    opts = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 30, 'Display', 'off', 'UseParallel', true);
else % 'pso'
    opts = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 30, 'Display', 'off', 'UseParallel', true);
end

for i=1:config.num_sample + 1
    
    ini_ctl_seq = repmat(cur_ctl(:,i), 1, config.T_sampling);
    
    if i==1
        current_sim_state = cur_state;
    else
        current_sim_state = state_sample_stack(:, i-1);
    end

    if strcmpi(solver_type, 'fmincon')
        % --- Your Original FMINCON Logic (Unchanged) ---
        ctl_seq_r = fmincon(@(ctl_seq) loss_fun(config, current_sim_state, ...
            [ini_ctl_seq(1:2*config.ctl_dim, :); ctl_seq], J2L_mat, "V_R",unit_ball_sample,...
            non_zero_count), ini_ctl_seq( 2*config.ctl_dim+1:end, :),[],[],...
            [],[],[],[], @(ctl_seq) mpc_cst(config, ctl_seq, "V_R"), opts);

        ctl_seq_t = fmincon(@(ctl_seq) loss_fun(config, current_sim_state, ...
            [ctl_seq; ctl_seq_r], J2L_mat, "V_T", unit_ball_sample, ...
            non_zero_count), ini_ctl_seq(1:2*config.ctl_dim, :), [],[],[], ...
            [],[],[], @(ctl_seq) mpc_cst(config, ctl_seq, "V_T"), opts);
        ctl_seq_t(config.ctl_dim + (1:config.ctl_dim), 1) = zeros([config.ctl_dim, 1]);

    else
        % --- New EA Solver Logic ---
        
        % --- RECON (Leader) SOLVER using EA ---
        recon_dims = config.num_rec_sat * config.ctl_dim * config.T_sampling;
        lb_r = -ones(recon_dims, 1); ub_r = ones(recon_dims, 1);
        % Use an anonymous function to pass extra arguments to the fitness function
        % obj_fun_r = @(x) ea_fitness_function(x, config, current_sim_state, ...
        %     ini_ctl_seq(1:2*config.ctl_dim, :), "V_R", J2L_mat, unit_ball_sample, non_zero_count);

        obj_fun_r = @(x) ea_fitness_function(x, config, current_sim_state, ...
            [zeros(config.ctl_dim,config.T_sampling); ini_ctl_seq(config.ctl_dim + (1:config.ctl_dim), :)], ...
            "V_R", J2L_mat, unit_ball_sample, non_zero_count);


        if strcmpi(solver_type, 'ga')
            ctl_seq_r_vec = ga(obj_fun_r, recon_dims, [],[],[],[],lb_r,ub_r,[],opts);
        else
            ctl_seq_r_vec = particleswarm(obj_fun_r, recon_dims, lb_r, ub_r, opts);
        end
        ctl_seq_r = reshape(ctl_seq_r_vec, [], config.T_sampling);

        % --- TARGET (Follower) SOLVER using EA ---
        target_dims = 2 *config.ctl_dim * config.T_sampling;
        lb_t = -ones(target_dims, 1); ub_t = ones(target_dims, 1);
        obj_fun_t = @(x) ea_fitness_function(x, config, current_sim_state, ...
            ctl_seq_r, "V_T", J2L_mat, unit_ball_sample, non_zero_count);

        if strcmpi(solver_type, 'ga')
            ctl_seq_t_vec = ga(obj_fun_t, target_dims, [],[],[],[],lb_t,ub_t,[],opts);
        else
            ctl_seq_t_vec = particleswarm(obj_fun_t, target_dims, lb_t, ub_t, opts);
        end
        ctl_seq_t = reshape(ctl_seq_t_vec, [], config.T_sampling);
    end
    
    ini_admissible_ctl(:, i) = [ctl_seq_t(:, 1); ctl_seq_r(:,1)];
    % ini_admissible_ctl(:, i) = [ctl_seq_t(:,1);ctl_seq_r(:,1)];

    if norm(ini_admissible_ctl(1:config.ctl_dim, i)) > config.tgt_att_u_max

            ini_admissible_ctl(1:config.ctl_dim, i) = config.tgt_att_u_max*...
            ini_admissible_ctl(1:config.ctl_dim, i)/norm(ini_admissible_ctl...
            (1:config.ctl_dim, i));
    end

    norm_u_max = config.tgt_orb_u_max;
    for k = 1:config.num_rec_sat+1
        if norm(ini_admissible_ctl(config.ctl_dim*k + (1:config.ctl_dim), i)) > ...
                config.tgt_orb_u_max
            if k > 1
                norm_u_max = config.rec_orb_u_max(k-1);
            end
            ini_admissible_ctl(config.ctl_dim*k + (1:config.ctl_dim), i) = ...
                norm_u_max*ini_admissible_ctl(config.ctl_dim*k +...
                (1:config.ctl_dim), i)/norm(ini_admissible_ctl(config.ctl_dim*k...
                + (1:config.ctl_dim), i));
        end
    end
end
% --- END OF REFACTORING ---

end