function [fun_value] = tgt_fun_po_ev(config, state, coeff)
%
fun_value = 0;
for i=1:size(state,2)
    fun_value = fun_value - coeff'*basis_fun(config.wgt_mat*state(:,i));
end

end