function [y] = smooth_operator(x,mu,u_bar,K)

if isfield(K,"l")
    num_l = K.l;
else
    num_l = 0;
end

t = u_bar - x/mu;

[d,P_cell] = eigK_with_vec(t,K);
n = length(d);
d_orig = d;

dsort = sort(d);
sumvalue = zeros(n,1) ;
lambda = 0 ;
for i = n-1: -1: 1
    lambda = dsort(i) ;
    g = (n-i)*(dsort(i+1) - dsort(i)) ;
    sumvalue(i) = sumvalue(i+1) + g ;
    if sumvalue(i) >=1
        lambda = (1-sumvalue(i))*((dsort(i+1)-dsort(i))/(sumvalue(i+1)-sumvalue(i)))+ dsort(i) ;
        break;
    end
    if (i == 1)
        lambda = dsort(1) + (sumvalue(1)-1)/n ;
    end
end



tmp=d_orig-lambda*ones(length(d),1);
tmp(tmp<0)=0;


[dim,~] = size(x);
y = zeros(dim,1);

if num_l > 0
    y(1:num_l) = tmp(1:num_l);
end

start = num_l; start_d = num_l;

for i = 1:length(K.s)
    tmp_n = K.s(i);
    Tmp = P_cell{i}*diag(tmp(start_d+1:start_d+tmp_n))*P_cell{i}.';
    y(start+1:start+tmp_n^2) = reshape(Tmp,[tmp_n^2,1]);

    start = start + tmp_n^2;
    start_d = start_d + tmp_n;
end


end

