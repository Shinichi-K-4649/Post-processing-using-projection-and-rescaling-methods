function [ppinfo,bound] = CheckInputSol(A,b,c,K,sol,skip)

x = sol.x; 
p_min = mineigK(x,K); 
ppinfo.x = []; bound.val = []; bound.sol = [];

if p_min > 0
    ppinfo.x = x; % Set the approximate optimal interior primal solution
else
    x_new = x + (abs(p_min) + 1e-15)*eyeK(K);
    if norm(A*x_new-b) <= 1e-4
        ppinfo.x = x; % Set the approximate optimal interior primal solution
    end
end

if skip == 0
    y = sol.y; z = sol.z;
    d_min1 = mineigK(z,K); d_min2 = mineigK(c-A.'*y,K);
    if d_min2 > 0
        bound.val = b'*y; bound.sol = y; % Set the current dual feasible solution
    else
        if d_min1 > 0
            tmp = -(z-c);
            [U,Sigma,V] = svd(full(A));
            tol = max(max(Sigma))/1e16;
            [m,~]=size(A);
            s = diag(Sigma(1:m,1:m));
            index = find(s<tol);
            de_rank = length(index);
            s = diag(Sigma(1:m-de_rank,1:m-de_rank)); s = 1./s;
            S = zeros(m,m-de_rank); S(1:m-de_rank,1:m-de_rank)=diag(s);
            y1 = U*S*V(:,1:m-de_rank)'*tmp;
            if mineigK(c-A'*y1,K) > 0
                bound.val = b'*y; bound.sol = y; % Set the current dual feasible solution
            end
        end
    end
end

end

