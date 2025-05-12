function [fun_value] = h(x)
%
if x < -0.01
    fun_value = -ones(size(x));
elseif -0.01 <= x & x <= 0.01 
    fun_value = (1-2*(x-0.01)/(0.02)).*(x+0.01).^2/(0.02)^2-...
    (1+2*(x+0.01)/(0.02)).*(x-0.01).^2/(0.02)^2;
else
    fun_value = ones(size(x));
end
end