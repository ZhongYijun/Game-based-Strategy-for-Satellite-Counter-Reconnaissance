function [coeff] = policy_evaluation(config, cur_state, cur_ctl, nxt_state,...
    coeff_ini, his_obs, J2L_mat, theta_max, trust_region_r, Hamilton_diff,...
    J_hamilton_fun)
%
coeff = zeros(size(coeff_ini));
opts = optimoptions('fmincon', 'Algorithm', 'sqp','MaxFunEvals', 5000);
%% 下面利用对偶性原理松弛求解值函数 V_T, V_R

coeff(:, 1) = fmincon(@(coeff) tgt_fun_po_ev(config, cur_state, coeff), ...
    coeff_ini(:, 1), [],[],[],[],[],[], @(coeff) quardic_cst(config, ...
    cur_state, cur_ctl, nxt_state, his_obs, J2L_mat, theta_max, coeff_ini(:,...
    1), coeff, trust_region_r, Hamilton_diff, J_hamilton_fun, "V_T"), opts);
coeff(:, 2) = fmincon(@(coeff) tgt_fun_po_ev(config, cur_state, coeff), ...
    coeff_ini(:, 2), [],[],[],[],[],[], @(coeff) quardic_cst(config, ...
    cur_state, cur_ctl, nxt_state, his_obs, J2L_mat, theta_max, coeff_ini(:,...
    2), coeff, trust_region_r, Hamilton_diff, J_hamilton_fun, "V_R"), opts);

end