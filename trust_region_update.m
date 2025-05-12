function [trust_region] = trust_region_update(config, trust_region,...
    real_hamilton_diff, pred_hamilton_diff)
%
p_k = norm(real_hamilton_diff)/norm(pred_hamilton_diff);
if p_k < 0.25
    trust_region = 0.25*trust_region;
elseif p_k > 0.75
    trust_region = min(2*trust_region, config.max_trust_redius);
end
end