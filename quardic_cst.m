function [ieq_cst,eq_cst] = quardic_cst(config, cur_state, cur_ctl, nxt_state,...
    his_obs, J2L_mat, theta_max, coeff_ini, coeff, trust_region_r, Hamilton_diff, ...
    J_hamilton_fun, type_value_fun)
%
ieq_cst = zeros([2*(config.num_sample +1) + 1, 1]);
eq_cst = [];
coeff_diff = coeff - coeff_ini;
for i=1:config.num_sample +1
    if type_value_fun == "V_T"

        ieq_cst(i,1) = Hamilton_fun(coeff_ini, cur_state(:, i), cur_ctl(:, i),...
            nxt_state(:, i), his_obs, J2L_mat, theta_max(:, i), config, "V_T") + ...
            coeff_diff'*Hamilton_diff{1, i} + 0.5*coeff_diff'*...
            J_hamilton_fun{1, i}*coeff_diff;
    else

        ieq_cst(i,1) = Hamilton_fun(coeff_ini, cur_state(:, i), cur_ctl(:, i),...
            nxt_state(:, i), his_obs, J2L_mat, theta_max(:, i), config, "V_R") + ...
            coeff_diff'*Hamilton_diff{2, i} + 0.5*coeff_diff'*...
            J_hamilton_fun{2, i}*coeff_diff;
    end
    ieq_cst(config.num_sample+ 1+ i, 1) = coeff_diff'*(basis_fun(config.wgt_mat...
        *cur_state(:, i)) - config.eta*basis_fun(config.wgt_mat*nxt_state(:, i)));
end 
if type_value_fun == "V_T"
    ieq_cst(end, 1) = trust_region_r(1) - norm(coeff_diff);
else
    ieq_cst(end, 1) = trust_region_r(2) - norm(coeff_diff);
end
ieq_cst = -ieq_cst;
end