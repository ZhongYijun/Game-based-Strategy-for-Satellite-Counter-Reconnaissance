function [his_obs] = his_obs_update(config, cur_r, his_obs)
%
%下面来更新历史观测区域面积
for i = 1:config.num_rec_sat
%若当前状态位于历史观测球外面时，视为有效观测区域增加，此时更新历史观策球
if dot(cur_r(:,i), his_obs(2:end, i)) < cos(his_obs(1, i))

    his_obs_ang_new = 0.5*(acos(dot(his_obs(2:end, i), cur_r(:,i)))+ ...
        his_obs(1, i));
    mu_ = dot(his_obs(2:end, i), cur_r(:,i));
    s_ = cos(his_obs_ang_new)/cos(his_obs_ang_new-his_obs(1, i));
    lambda_1 = (1-mu_*s_)/(1+s_)/(1-mu_);
    lambda_2 = 1-lambda_1;
    his_obs(2:end, i) = (lambda_1*his_obs(2:end,i) + lambda_2*cur_r(:,i))/...
        norm(lambda_1*his_obs(2:end,i) + lambda_2*cur_r(:,i));

end

end

end