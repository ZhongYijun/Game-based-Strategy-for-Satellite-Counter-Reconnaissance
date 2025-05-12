function [trans_mat, cur_r] = att_trans_mat(config, cur_state)
%
sigma_=[0, -cur_state(3,1), cur_state(2,1); cur_state(3,1), 0, -cur_state(1,1);...
    -cur_state(2,1), cur_state(1,1), 0];
trans_mat = eye(3)-4*(1-norm(cur_state(1:3,1))^2)/(1+norm(cur_state(1:3,1))...
    ^2)^2*sigma_+8/(1+norm(cur_state(1:3,1))^2)^2*sigma_*sigma_;
cur_r = zeros(3, config.num_rec_sat);
for i=1:config.num_rec_sat
    cur_r_i = cur_state((i-1)*config.orb_dim+config.att_dim+(1:3), 1);
    cur_r(:, i) = cur_r_i/norm(cur_r_i);
end
end