function [A_aft,K_aft] = transform_type1(A,b,c,K,InputValue)
% Construct the primal model

[m,~] = size(A);
if isfield(K,'l') == 1
    num_l = K.l;
else
    num_l = 0;
end

A_aft = [-b,zeros(m,1),A]; tmp = [-InputValue,1,c.'];
A_aft = [A_aft;tmp];

K_aft = K; K_aft.l = num_l + 2;

end

