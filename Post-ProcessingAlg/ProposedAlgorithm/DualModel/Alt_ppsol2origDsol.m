function [solinfo] = Alt_ppsol2origDsol(solinfo,A,c)
% Construct a dual solution using the outputs from the projection and rescaling method

num = length(solinfo.sol);
kappa = solinfo.sol(num); y = solinfo.sol(1:num-1);

if kappa < 0
    orig_sol = []; solinfo.type="unstable output dual";
else
    orig_sol.y = -y/kappa;
    orig_sol.z = c - A.'*orig_sol.y;
end

solinfo.origsol = orig_sol;

end

