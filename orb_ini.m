function [eqn] = orb_ini(w_ref,cur_state)
% 函数用于生成绕飞轨道初始初始状态
eqn=zeros([2,1]);
eqn(1,1) = sqrt((2*cur_state(6,1)/w_ref-3*cur_state(1,1))^2+(cur_state(4,1)...
            /w_ref)^2)-2500;
eqn(2,1) = cur_state(6,1)-2*w_ref*cur_state(1,1);
end