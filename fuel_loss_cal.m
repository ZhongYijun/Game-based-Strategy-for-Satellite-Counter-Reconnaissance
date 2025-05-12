function [eng_loss_T, eng_loss_R] = fuel_loss_cal(cur_ctl_att_tgt,...
    cur_ctl_orb_tgt, cur_ctl_orb_rec)
%函数用于计算侦察卫星与目标卫星的燃料消耗
eng_loss_T = zeros(2,1);
eng_loss_T(1,1) = norm(cur_ctl_att_tgt);
eng_loss_T(2,1) = norm(cur_ctl_orb_tgt);
eng_loss_R = norm(cur_ctl_orb_rec);
end