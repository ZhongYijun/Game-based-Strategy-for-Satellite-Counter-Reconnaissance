% Filename: ea_fitness_function.m

function fitness = ea_fitness_function(ctl_seq_vector, config, state, ...
    fixed_opponent_ctl_seq, type_value_fun, J2L_mat, unit_ball_sample, non_zero_count)
%
% This is the ADAPTER or "Fitness Function" for MATLAB's GA and PSO solvers.
% It takes a 1D vector of control variables, reshapes it, and returns a single
% scalar cost by combining the original loss with a penalty for constraint violations.
%

%% 1. Reshape the 1D input vector back to a 2D control sequence matrix
ctl_seq_optim = reshape(ctl_seq_vector, [], config.T_sampling);

%% 2. Assemble the full control sequence for simulation
if strcmp(type_value_fun, "V_T")
    ctl_seq_optim(1:config.ctl_dim, :) = zeros(config.ctl_dim, config.T_sampling);
    % We are optimizing Target's control, Recon's is fixed
    full_ctl_seq = [ctl_seq_optim; fixed_opponent_ctl_seq];
else % "V_R"
    % We are optimizing Recon's control, Target's is fixed
    full_ctl_seq = [fixed_opponent_ctl_seq; ctl_seq_optim];
end

%% 3. Calculate the primary loss from your original function
primary_loss = loss_fun(config, state, full_ctl_seq, J2L_mat, ...
    type_value_fun, unit_ball_sample, non_zero_count);

%% 4. Calculate the penalty for constraint violations
% Get constraint values from your original constraint function
[ineq_cst, eq_cst] = mpc_cst(config, ctl_seq_optim, type_value_fun);

% A large weight to heavily penalize any constraint violation
penalty_weight = 1e8; 

% For inequality constraints (g(x) <= 0), the penalty is for any g(x) > 0.
penalty = sum(max(0, ineq_cst).^2);

% For equality constraints (h(x) == 0), the penalty is for any h(x) != 0.
if ~isempty(eq_cst)
    penalty = penalty + sum(eq_cst.^2);
end

%% 5. Combine into a single fitness value for the solver
fitness = primary_loss + penalty_weight * penalty;

end