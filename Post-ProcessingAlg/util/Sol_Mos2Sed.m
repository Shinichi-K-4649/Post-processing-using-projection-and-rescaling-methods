function [psol,dsol] = Sol_Mos2Sed(prob,res)

psol=[];start =0;
Psd = prob.bardim;
dsol.y = res.sol.itr.y;
dsol.z = [];
for i = 1:length(Psd)
    n_s = Psd(i);
    X = ones(n_s,n_s); X = tril(X);
    I = find(X);
    num = (n_s+1)*n_s/2;

    X(I) = res.sol.itr.barx(start+1:start+num);
    X = X+X.' - diag(diag(X));
    psol = [psol;vec(X)];

    X = sparse(n_s,n_s);
    X(I) = res.sol.itr.bars(start+1:start+num);
    X = X+X.' - diag(diag(X));
    dsol.z = [dsol.z;vec(X)];

    start = start + num;
end


if isfield(res.sol.itr,"xx")
    psol_l = res.sol.itr.xx;
    psol = [psol_l;psol];

    dsol_l = res.sol.itr.slx;
    dsol.z = [dsol_l;dsol.z];
end



end

