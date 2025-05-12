function [co_obser_area] = co_obser_area_cal(config, his_obs, J2L_mat,...
    theta_max)
%
% 输入：当前状态: cur_state 历史观测区域近似: his_obs

co_obser_area = integral2(@(psi,theta) obser_fun(config, psi, theta,...
    his_obs, J2L_mat, theta_max), 0, pi,0, 2*pi,'AbsTol', 1e-6, 'RelTol', 1e-6);


end