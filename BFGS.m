function [H_k] = BFGS(y_k, s_k, H_k)
%
L_k = (eye(size(y_k,1)) - s_k*y_k'/(s_k'*y_k));
H_k = L_k*H_k*L_k'+ s_k*s_k'/(s_k'*y_k);
end