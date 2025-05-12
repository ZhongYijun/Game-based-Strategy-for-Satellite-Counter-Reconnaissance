function [ctl_att_tgt] = u_costate_att(config, costate)
%
% 
 if norm(costate(4:config.att_dim))> 2*config.tgt_att_u_max*config.wgt_T2

        ctl_att_tgt = -config.tgt_att_u_max*costate(4:config.att_dim)/...
            norm(costate(4:config.att_dim));

 else

        ctl_att_tgt = -0.5*costate(4:config.att_dim)/config.wgt_T2;
 end
end