function [nxt_state] = Simulator(config, cur_state, cur_ctl)
%% Dynamic model of attitude and orbit
%  Input:
%  Output:
cur_ctl_att_tgt = cur_ctl(1:config.ctl_dim, 1);
cur_ctl_orb_tgt = cur_ctl(config.ctl_dim+1:2*config.ctl_dim, 1);
cur_ctl_orb_rec = cur_ctl(2*config.ctl_dim+1: end, 1);
nxt_state = zeros(size(cur_state));
sigma_ = [0, -cur_state(3),cur_state(2);cur_state(3),0,-cur_state(1);
    -cur_state(2),cur_state(1),0];
omega_ = [0, -cur_state(6),cur_state(5);cur_state(6),0,-cur_state(4);
    -cur_state(5),cur_state(4),0];
M_sigma = 0.25*((1-norm(cur_state(1:3))^2)*eye(3)+2*sigma_+2*...
    cur_state(1:3)*cur_state(1:3)');
nxt_state(1:3) = cur_state(1:3)+...
    config.sample_time_span*M_sigma*cur_state(4:6);
nxt_state(4:6) = cur_state(4:6)+config.sample_time_span*config.MI\(-omega_*config.MI...
    *cur_state(4:6)+cur_ctl_att_tgt);
    nxt_state(config.att_dim+1:end) = kron...
        (eye(config.num_rec_sat),config.A_ZOH)*cur_state(config.att_dim...
        +1: end) + kron(eye(config.num_rec_sat),config.B_ZOH)*...
        (cur_ctl_orb_rec-kron(ones([config.num_rec_sat,1]),cur_ctl_orb_tgt));
sigma = norm(nxt_state(1:3, 1));
if sigma>1
    nxt_state(1:3, 1) = -nxt_state(1:3, 1)/sigma^2;
end
end