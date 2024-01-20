function [solinfo] = Alt_ppsol2origPsol(solinfo,c,K)
% Construct a primal solution, a reducing direction or an improving ray using the outputs from the projection and rescaling method

[d,~] = size(solinfo.sol);
tau = solinfo.sol(1); x = solinfo.sol(3:d);

if abs(tau) <= 1e-12
    value1 = c'*x;
    if abs(value1) <= 1e-12
        value2 = mineigK(x,K);
        value3 = norm(x);
        if value3 > 1e-12 && value2 >= -1e-12
            solinfo.type = "dual reducing direction";
            orig_sol = x;
        end
    end

    if value1 < -1e-12
        tmp = -x/value1;
        value2 = mineigK(tmp,K);
        if value2 >= -1e-12 
            solinfo.type = "primal improving ray";
            orig_sol = tmp;
        end
    end
end


if solinfo.type ~= "dual reducing direction" && solinfo.type ~= "primal improving ray"
    if tau < 0
        orig_sol = []; solinfo.type="unstable output primal";
    else
        orig_sol = x/tau;
    end
end

solinfo.origsol = orig_sol;

end

