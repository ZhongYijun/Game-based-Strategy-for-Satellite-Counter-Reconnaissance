function [ctl_orb_rec] = u_costate_orb_rec(config,costate)
%
% Input: costate
 if norm(config.B_ZOH'*costate)> 2*config.sample_time_span*config.rec_orb_u_max*...
         config.wgt_R2

        ctl_orb_rec = config.rec_orb_u_max*config.B_ZOH'*costate/...
            norm(config.B_ZOH'*costate);

 else

        ctl_orb_rec = 0.5*config.B_ZOH'*costate/config.sample_time_span/...
            config.wgt_R2;
 end
end