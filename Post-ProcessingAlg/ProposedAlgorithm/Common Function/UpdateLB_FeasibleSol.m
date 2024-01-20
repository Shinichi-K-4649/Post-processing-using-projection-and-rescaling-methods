function [bound] = UpdateLB_FeasibleSol(A,b,c,K,bound,y2)
% Function to create a better dual feasible solution using two dual feasible solutions

y1 = bound.sol; update_flag = 1; delta = b.'*(y1-y2);
if delta > 0
    d = y1-y2; orig = y1;
elseif delta < 0
    d = y2-y1; orig = y2; delta = -delta; bound.sol = y2; bound.val = b.'*y2;
else
    update_flag = 0;
end

if update_flag == 1
    low = 0; alpha = 5;
    while true
        sol = orig + alpha*d;
        Slack = c-A.'*sol;
        if mineigK(Slack,K) > 0
            bound.sol = sol; bound.val = b.'*sol;
            break
        else
            up = alpha;
            alpha = (up+low)/2;
        end

        if abs(alpha*delta) <= 1e-16 || alpha <= 1e-12
            break
        end
    end

end

end

