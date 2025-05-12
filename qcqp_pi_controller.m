function [cur_ctl] = qcqp_pi_controller(config, cur_state, state_sample_stack,...
    his_obs, cur_ctl, J2L_mat)
% 基于QCQP约束优化方法的PI迭代算法


lst_coeff_vec = zeros([1+size(cur_state, 1)*(size(cur_state, 1)+3)/2, 2]);
trust_region_r = 0.5*config.max_trust_redius*ones([2,1]);
% 初始化Hamilton函数对应的Jacobian矩阵与梯度向量
J_hamilton_fun = cell([2, config.num_sample + 1]);
cur_state = [cur_state, state_sample_stack];
Hamilton_diff = cell([2, config.num_sample + 1]);
pred_hamilton_diff = [0;0];
real_hamilton_diff = [0;0];
for j=1:config.num_sample + 1
    J_hamilton_fun{1, j} = eye(size(lst_coeff_vec, 1));
    J_hamilton_fun{2, j} = eye(size(lst_coeff_vec, 1));
end

theta_max = zeros([config.num_rec_sat, config.num_sample + 1]);
for i=1:config.num_rec_sat
for j=1:config.num_sample + 1
    theta_max(i, j) = acos(config.radius_tgt_sat/norm(cur_state((i-1)*...
          config.orb_dim+config.att_dim+(1:3), j)));
end
nxt_state = zeros(size(cur_state));

for i=1:config.max_policy_iter
    for j=1:config.num_sample + 1
        nxt_state(:, j) = Simulator(config, cur_state(:, j), cur_ctl(:, j));
        [Hamilton_diff{1, j}, Hamilton_diff{2, j}] = hamilton_diff(cur_state...
            (:, j), cur_ctl(:, j), nxt_state(:, j), config, lst_coeff_vec);

    end
    
    coeff_vec = policy_evaluation(config, cur_state, cur_ctl, nxt_state,...
        lst_coeff_vec, his_obs, J2L_mat, theta_max, trust_region_r, Hamilton_diff,...
        J_hamilton_fun);

    cur_ctl = policy_improvement(config, cur_state, his_obs, J2L_mat, ...
        theta_max, coeff_vec);

    coeff_diff = coeff_vec - lst_coeff_vec;

    % 对于Hamilton函数对应的Jacobian矩阵进行修正
    % for j=1:size(cur_state,1)
    %     grad_coeff_ini = gradient(@(coeff) Hamilton_fun(coeff, cur_state(:,...
    %         j), nxt_state(:, i), his_obs, config), coeff_ini);
    %     grad_coeff = gradient(@(coeff) Hamilton_fun(coeff, cur_state(:, j),...
    %     nxt_state(:, i), his_obs, J2L_mat, config), coeff_ini);
    %     Hamilon_diff = grad_coeff - grad_coeff_ini;
    for j=1:config.num_sample + 1
        pred_hamilton_diff(1) = pred_hamilton_diff(1) + coeff_diff(:,1)'*...
            Hamilton_diff{1, j} + 0.5*coeff_diff(:,1)'*J_hamilton_fun{1, j}...
            *coeff_diff(:,1);

        pred_hamilton_diff(2) = pred_hamilton_diff(2) + coeff_diff(:,2)'*...
            Hamilton_diff{2, j} + 0.5*coeff_diff(:,2)'*J_hamilton_fun{2, j}...
            *coeff_diff(:,2);

        real_hamilton_diff(1) = real_hamilton_diff(1) + Hamilton_fun(coeff_vec,...
            cur_state(:, j), cur_ctl(:, j), nxt_state(:, j), his_obs, J2L_mat,...
            theta_max(:, j), config, "V_T")-Hamilton_fun(lst_coeff_vec,...
            cur_state(:, j), cur_ctl(:, j), nxt_state(:, j), his_obs, J2L_mat,...
            theta_max(:, j), config, "V_T"); 

        real_hamilton_diff(2) = real_hamilton_diff(2) + Hamilton_fun(coeff_vec,...
            cur_state(:, j), cur_ctl(:, j), nxt_state(:, j), his_obs, J2L_mat,...
            theta_max(:, j), config, "V_R")-Hamilton_fun(lst_coeff_vec,...
            cur_state(:, j), cur_ctl(:, j), nxt_state(:, j), his_obs, J2L_mat,...
            theta_max(:, j), config, "V_R"); 

        J_hamilton_fun{1, j} = BFGS(coeff_diff(:, 1), Hamilton_diff{1, j},...
            J_hamilton_fun{1, j});
        J_hamilton_fun{2, j} = BFGS(coeff_diff(:, 2), Hamilton_diff{2, j},...
            J_hamilton_fun{2, j});
    end
    % 
    trust_region_r(1) = trust_region_update(config, trust_region_r(1), ...
    real_hamilton_diff(1), pred_hamilton_diff(1));

    trust_region_r(2) = trust_region_update(config, trust_region_r(2), ...
    real_hamilton_diff(2), pred_hamilton_diff(2));

    if norm(coeff_diff) < config.err_tol
        break
    end
    lst_coeff_vec = coeff_vec;
    % grad_coeff_ini = ;
end
end