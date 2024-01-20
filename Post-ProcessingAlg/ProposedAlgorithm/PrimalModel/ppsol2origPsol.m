function [solinfo] = ppsol2origPsol(solinfo)
% Construct a primal solution using the outputs from the projection and rescaling method
[d,~] = size(solinfo.sol);
tau = solinfo.sol(1); x = solinfo.sol(3:d);

if tau < 0
    solinfo.type = "unstable output primal";
end

orig_sol = x/tau;
solinfo.origsol = orig_sol;
end

