function [loc] = unit_ball_sample_gen(config)
%

nums_points = config.max_num_ball_sample; % 点数量
radius = 1;
loc = zeros(nums_points, 3);
for i = 1 : nums_points
    phi = acos(-1.0 + (2.0 * i - 1.0) / nums_points);
    theta = sqrt(nums_points * pi) * phi;
    loc(i, 1) = radius * cos(theta) * sin(phi);
    loc(i, 2) = radius * sin(theta) * sin(phi);
    loc(i, 3) = radius * cos(phi);
end
end