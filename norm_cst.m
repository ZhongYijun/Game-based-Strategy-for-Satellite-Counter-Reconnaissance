function [ieq_cst, eq_cst] = norm_cst(config, cur_state, cur_r, tgt_r, J2L_mat, trans_mat,...
    type_value_fun)
%
if size(tgt_r,1)~=3
    tgt_r = tgt_r';
end
eq_cst = zeros([2*config.num_rec_sat, 1]);
ieq_cst = zeros([2*config.num_rec_sat, 1]);
max_norm =cos(pi/8)*ones([config.num_rec_sat, 1]) ;
% for i=0:config.T_sampling-1
%     for j=1:config.num_rec_sat
%         cur_rel_orb = (config.A_ZOH^i)*config.B_ZOH;
%         max_norm(j) = max_norm(j)  - config.rec_orb_u_max(j)*norm(cur_r(:,j)'...
%             *trans_mat*cur_rel_orb(1:3,:));
%     end
% end
% 
% for j=1:config.num_rec_sat
%     cur_rel_orb = (config.A_ZOH^config.T_sampling)*cur_state(config.att_dim...
%         +(j-1)*config.orb_dim + (1:config.orb_ctl_dim), 1);
%     max_norm(j) = max_norm(j) + cur_r(:,j)'*trans_mat*cur_rel_orb;
% end

if type_value_fun == "V_R"
    for i=1:config.num_rec_sat
        ieq_cst(i,1) = tgt_r(:, i)'*J2L_mat*config.a_sun-cos(config.ang_sun);
        % if cur_r'*his_obs(2:end, i) < his_obs(1, i)
        %     ieq_cst(2*i, 1) = min(max_norm(i), his_obs) - cur_r(:, i)'*...
        %         tgt_r(:, i);
        % else
        if cur_r(:, i)'*J2L_mat*config.a_sun > cos(config.ang_sun)
            eq_cst(i+config.num_rec_sat, 1) = tgt_r(:, i)'*J2L_mat*...
                config.a_sun - cos(config.ang_sun); 
        else
            ieq_cst(i+config.num_rec_sat, 1) = max_norm(i) - cur_r(:, i)'*...
                tgt_r(:, i);
        end
        % end
        eq_cst(i,1) = norm(tgt_r(:, i)) - 1;
    
    end
else
    for i=1:config.num_rec_sat
        ieq_cst(i, 1) = max_norm(i) - cur_r(:, i)'*tgt_r(:, i);
        eq_cst(i,1) = norm(tgt_r(:, i)) - 1;
    end

end

end