function [y] = h2(x,lower_bnd)
%
if x<lower_bnd
    y = 1;
elseif lower_bnd <= x && x <= lower_bnd + 0.01 
    y = (1+2*(x-lower_bnd)/(0.01)).*(x-lower_bnd- 0.01).^2/(0.01)^2;
else
    y=0;
end
end