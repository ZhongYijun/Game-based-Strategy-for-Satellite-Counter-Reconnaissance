function [ieq_cst, eq_cst] = po_im_cst(config, cur_ctl, type_value_fun)
%
eq_cst = [];
if type_value_fun == "V_T"
    ieq_cst = zeros([2, 1]);
    ieq_cst(1, 1) = norm(cur_ctl(1:config.ctl_dim))-config.tgt_att_u_max;
    ieq_cst(2, 1) = norm(cur_ctl(config.ctl_dim+(1:config.ctl_dim)))-...
            config.tgt_orb_u_max;
else
    ieq_cst = zeros([config.num_rec_sat, 1]);
    for i=1:config.num_rec_sat
        ieq_cst(i, 1) = norm(cur_ctl((i-1)*config.ctl_dim+(1:config.ctl_dim)))...
             - config.rec_orb_u_max(i);
    end
end
end