function [solinfo] = ppsol2origDsol(solinfo,A,c,b,K)
% Construct a dual solution, a reducing direction or an improving ray using the outputs from the projection and rescaling method

num = length(solinfo.sol);
kappa = solinfo.sol(num); y = solinfo.sol(1:num-1);

if abs(kappa) <= 1e-12
    value1 = -b.'*y;
    z = A'*y; value2 = mineigK(z,K);
    if abs(value1) <= 1e-12 && value2 > -1e-12
        value3 = norm(z);
        if value3 > 1e-12
            solinfo.type = "primal reducing direction";
            orig_sol.y = y;
        end
    end

    if value1 > 1e-12
        tmp = y/value1;
        z = A'*tmp; value2 = mineigK(z,K);
        if value2 >= -1e-12
            solinfo.type = "dual improving ray";
            orig_sol.y = tmp;
        end
    end
end

if solinfo.type ~= "primal reducing direction" && solinfo.type ~= "dual improving ray"
    if kappa < 0
        orig_sol = []; solinfo.type="unstable output dual";
    else
        orig_sol.y = -y/kappa;
        orig_sol.z = c - A.'*orig_sol.y;
    end
end

solinfo.origsol = orig_sol;

end

