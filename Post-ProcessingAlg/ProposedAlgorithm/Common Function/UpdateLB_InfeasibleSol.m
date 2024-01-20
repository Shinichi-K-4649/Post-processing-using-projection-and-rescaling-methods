function [bound] = UpdateLB_InfeasibleSol(A,b,c,K,bound,y)
% Function to create a better dual feasible solution using dual feasible and infeasible solutions

update_flag = 0;
if b.'*y >= bound.val
    alpha = 2;
    d = y - bound.sol;
    while true
        new_sol = bound.sol + d/alpha;
        if mineigK(c - A.'*new_sol,K) > 0
            update_flag=1;
        end

        if update_flag == 1 && 1/alpha <= 1e-6 
            up =1/sqrt(alpha); 
            low = 1/alpha; 
            while true
                beta = up*0.9 + low*0.1;
                new_sol = bound.sol + d*beta;
                if mineigK(c - A.'*new_sol,K) > 0
                    break
                end
                up = up*0.1;
                if up <= low
                    break
                end
            end
        end

        if update_flag == 1 
            bound.sol = new_sol; bound.val = b.'*new_sol;
            break
        end

        alpha = alpha*alpha; 
        if abs(b'*d/alpha) <= 1e-16 || 1/alpha <= 1e-12
            break
        end
    end

else 
    d = bound.sol - y; 
    delta = b.'*d;
    up = 1; low = 0; 
    alpha = up;
    while true
        sol = bound.sol + alpha*d;
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

