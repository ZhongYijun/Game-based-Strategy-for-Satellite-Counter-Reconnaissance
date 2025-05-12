function [basis_f] = basis_fun(cur_state)
%

basis_f = [1;cur_state;];
for i=1:size(cur_state,1)
    basis_f = [basis_f; 0.5*cur_state(i)^2];
    for j=i+1:size(cur_state,1)
        basis_f = [basis_f; cur_state(i)*cur_state(j)];
    end
end
end