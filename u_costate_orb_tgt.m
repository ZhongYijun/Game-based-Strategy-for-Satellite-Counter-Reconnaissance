function [ctl_orb_tgt] = u_costate_orb_tgt(config, costate)
%
costate_ = zeros([config.orb_dim, 1]);
for i=1: config.num_rec_sat
costate_ = costate_ + costate((i-1)*config.orb_dim + 1:i*config.orb_dim, 1);
 if norm(costate_)> 2*config.sample_time_span*config.tgt_orb_u_max*...
         config.wgt_T2

        ctl_orb_tgt = config.tgt_orb_u_max*config.B_ZOH'*costate_/...
            norm(config.B_ZOH'*costate_);

 else

        ctl_orb_tgt = 0.5*config.B_ZOH'*costate_/config.sample_time_span...
            /config.wgt_T2;
 end
end