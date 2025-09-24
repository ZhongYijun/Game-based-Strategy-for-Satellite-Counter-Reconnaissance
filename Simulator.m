% function [nxt_state] = Simulator(config, cur_state, cur_ctl)
% %% Dynamic model of attitude and orbit
% %  Input:
% %  Output:
% cur_ctl_att_tgt = cur_ctl(1:config.ctl_dim, 1);
% cur_ctl_orb_tgt = cur_ctl(config.ctl_dim+1:2*config.ctl_dim, 1);
% cur_ctl_orb_rec = cur_ctl(2*config.ctl_dim+1: end, 1);
% nxt_state = zeros(size(cur_state));
% sigma_ = [0, -cur_state(3),cur_state(2);cur_state(3),0,-cur_state(1);
%     -cur_state(2),cur_state(1),0];
% omega_ = [0, -cur_state(6),cur_state(5);cur_state(6),0,-cur_state(4);
%     -cur_state(5),cur_state(4),0];
% M_sigma = 0.25*((1-norm(cur_state(1:3))^2)*eye(3)+2*sigma_+2*...
%     cur_state(1:3)*cur_state(1:3)');
% nxt_state(1:3) = cur_state(1:3)+...
%     config.sample_time_span*M_sigma*cur_state(4:6);
% nxt_state(4:6) = cur_state(4:6)+config.sample_time_span*config.MI\(-omega_*config.MI...
%     *cur_state(4:6)+cur_ctl_att_tgt);
%     nxt_state(config.att_dim+1:end) = kron...
%         (eye(config.num_rec_sat),config.A_ZOH)*cur_state(config.att_dim...
%         +1: end) + kron(eye(config.num_rec_sat),config.B_ZOH)*...
%         (cur_ctl_orb_rec-kron(ones([config.num_rec_sat,1]),cur_ctl_orb_tgt));
% sigma = norm(nxt_state(1:3, 1));
% if sigma>1
%     nxt_state(1:3, 1) = -nxt_state(1:3, 1)/sigma^2;
% end
% end

% Filename: Simulator.m (CORRECTED & FINAL VERSION)

function [nxt_state] = Simulator(config, cur_state, cur_ctl)
%% Dynamic model of attitude and orbit (with corrected attitude kinematics)

cur_ctl_att_tgt = cur_ctl(1:config.ctl_dim, 1);
cur_ctl_orb_tgt = cur_ctl(config.ctl_dim+1:2*config.ctl_dim, 1);
% --- FIX for 1-on-1 scenario ---
cur_ctl_orb_rec = cur_ctl(2*config.ctl_dim+1: end, 1);
if size(cur_ctl_orb_rec, 1) > config.ctl_dim * config.num_rec_sat
    % This handles the case where the input vector might be larger
    cur_ctl_orb_rec = cur_ctl_orb_rec(1 : config.ctl_dim * config.num_rec_sat);
end
% --- END FIX ---

nxt_state = zeros(size(cur_state));

%% Attitude Dynamics (Inertial angular velocity propagation) - This part is correct
omega_vec = cur_state(4:6);
omega_mat = [0, -omega_vec(3), omega_vec(2); omega_vec(3), 0, -omega_vec(1); -omega_vec(2), omega_vec(1), 0];
nxt_state(4:6) = omega_vec + config.sample_time_span * (config.MI \ (-omega_mat * config.MI * omega_vec + cur_ctl_att_tgt));

%% Attitude Kinematics (MRP propagation) - This part requires correction
sigma_vec = cur_state(1:3);
sigma_mat = [0, -sigma_vec(3), sigma_vec(2); sigma_vec(3), 0, -sigma_vec(1); -sigma_vec(2), sigma_vec(1), 0];
M_sigma = 0.25 * ((1 - sigma_vec'*sigma_vec) * eye(3) + 2*sigma_mat + 2 * (sigma_vec * sigma_vec'));

% --- START OF BUG FIX ---
% To propagate the LVLH-to-Body attitude, we need the RELATIVE angular velocity.
% omega_rel = omega_inertial - C_body_LVLH * omega_LVLH
% where C_body_LVLH is the transpose of trans_mat (LVLH_to_Body).

% 1. Calculate the current LVLH-to-Body transformation matrix
trans_mat = att_trans_mat_local(sigma_vec);

% 2. Define LVLH angular velocity in LVLH frame (rotation is around Z-axis)
omega_LVLH_in_LVLH = [0; 0; config.orb_w];

% 3. Calculate the relative angular velocity, expressed in the BODY frame
omega_inertial_in_body = cur_state(4:6);
omega_rel_in_body = omega_inertial_in_body - trans_mat * omega_LVLH_in_LVLH;

% 4. Use the CORRECT relative angular velocity to propagate MRPs
nxt_state(1:3) = sigma_vec + config.sample_time_span * M_sigma * omega_rel_in_body;
% --- END OF BUG FIX ---

%% Orbit Dynamics (Unchanged)
nxt_state(config.att_dim+1:end) = kron(eye(config.num_rec_sat), config.A_ZOH) * cur_state(config.att_dim+1:end) ...
    + kron(eye(config.num_rec_sat), config.B_ZOH) * (cur_ctl_orb_rec - kron(ones([config.num_rec_sat,1]), cur_ctl_orb_tgt));

%% MRP Shadow Set Switching (Unchanged)
sigma_norm = norm(nxt_state(1:3, 1));
if sigma_norm > 1
    nxt_state(1:3, 1) = -nxt_state(1:3, 1) / (sigma_norm^2);
end

end

% --- Helper function to avoid dependency issues ---
function trans_mat = att_trans_mat_local(sigma_vec)
    sigma_norm_sq = sigma_vec' * sigma_vec;
    sigma_mat = [0, -sigma_vec(3), sigma_vec(2); sigma_vec(3), 0, -sigma_vec(1); -sigma_vec(2), sigma_vec(1), 0];
    trans_mat = eye(3) + ( -4*(1-sigma_norm_sq)*sigma_mat + 8*sigma_mat*sigma_mat ) / (1+sigma_norm_sq)^2;
end