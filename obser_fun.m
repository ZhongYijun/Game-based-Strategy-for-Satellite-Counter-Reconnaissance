function [int_fun] = obser_fun(config, psi, theta, his_obs, J2L_mat, theta_max)
%% 函数用于计算当前时刻的有效观测面积（考虑光照影响）
x = sin(psi).*cos(theta);
y = sin(psi).*sin(theta);
z = cos(psi);
his_obs = his_obs_update(his_obs);
% LVLH轨道坐标系下太阳方向向量
a_sun_OC = J2L_mat*config.a_sun;
int_fun = ones(size(psi));
int_fun = 0.5*int_fun.*(1+h(cos(min(theta_max(1), 0.5*(pi-...
        config.camera_view)))- x*his_obs(2,1) - y*his_obs(3,1) - z*his_obs(4,1)));
int_fun = 0.5*int_fun.*(1+h(x*a_sun_OC(1) + y*a_sun_OC(2) + z*a_sun_OC(3)...
    -cos(config.ang_sun)));

%% 利用分段样条插值方法，对于函数奇异性进行处理以避免积分不收敛问题
for i = 2:config.num_rec_sat
    int_fun_i = ones(size(psi));
    int_fun_i = 0.5*int_fun_i.*(1+h(cos(min(theta_max(i), 0.5*(pi-...
        config.camera_view)))- x*his_obs(2,i) - y*his_obs(3,i) - z*his_obs(4,i)));
    int_fun_i = 0.5*int_fun_i.*(1+h(x*a_sun_OC(1) + y*a_sun_OC(2) + z*a_sun_OC(3)...
    -cos(config.ang_sun)));
    int_fun = max(int_fun, int_fun_i);
end
int_fun = int_fun.*sin(psi);
end