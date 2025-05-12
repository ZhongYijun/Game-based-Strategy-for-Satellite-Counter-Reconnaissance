function [integral_val] = tgt_state_generator(config, psi, theta, his_obs,...
            cur_r, J2L_mat, theta_max, tgt_r)
%
x = sin(psi).*cos(theta);
y = sin(psi).*sin(theta);
z = cos(psi);
% his_obs = his_obs_update(config, tgt_r, his_obs);

int_fun = zeros(size(psi));

%% 首先来计算当前状态对应的观测范围
for i = 1:config.num_rec_sat
    int_fun_i = ones(size(psi));
    int_fun_i = 0.5*int_fun_i.*(1+h(cos(min(theta_max(1), (0.5*pi-...
        config.camera_view)))- x*tgt_r(1, i) - y*tgt_r(2, i) - ...
        z*tgt_r(3, i)));
    int_fun = max(int_fun, int_fun_i);
end

%% 接下来计算历史状态对应的观测范围
for i = 1:config.num_rec_sat

    int_fun_i = ones(size(psi));
    int_fun_i = 0.5*int_fun_i.*(1+h(cos(his_obs(1,i))- x*his_obs(2,i) - y*his_obs(3,i) - ...
        z*his_obs(4,i)));

    int_fun = max(int_fun, int_fun_i);
end
%% 最后来计算考虑光照条件下的观测范围
for i = 1:config.num_rec_sat
    a_sun = J2L_mat*config.a_sun;
    int_fun = 0.5*int_fun.*(1+h(x*a_sun(1) + y*a_sun(2) + z*a_sun(3)...
        -cos(config.ang_sun)));
end

integral_val = int_fun.*sin(psi);
end