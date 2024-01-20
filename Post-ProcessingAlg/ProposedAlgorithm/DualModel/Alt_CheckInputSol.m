function [ppinfo,bound,StockD] = Alt_CheckInputSol(A,b,c,K,sol)

y = sol.y; z = sol.z;
d_min1 = mineigK(z,K); d_min2 = mineigK(c-A.'*y,K);

ppinfo.z = []; bound.val = []; bound.sol = [];
StockD.Sol = []; StockD.err4 = [];

if d_min1 > 0
    ppinfo.z = z; % Set the approximate optimal interior dual solution
else
    z_new = z + (abs(d_min1) + 1e-15)*eyeK(K);
    if norm(c-A'*y-z_new) <= 1e-4
        ppinfo.z = z; % Set the approximate optimal interior dual solution
    end
end

if ~isempty(ppinfo.z)
    tmp = -(ppinfo.z-c);
    [U,Sigma,V] = svd(full(A));
    tol = max(max(Sigma))/1e16;
    [m,~]=size(A);
    s = diag(Sigma(1:m,1:m));
    index = find(s<tol);
    de_rank = length(index);
    s = diag(Sigma(1:m-de_rank,1:m-de_rank)); s = 1./s;
    S = zeros(m,m-de_rank); S(1:m-de_rank,1:m-de_rank)=diag(s);
    y1 = U*S*V(:,1:m-de_rank)'*tmp;

    ppinfo.y = y1;
    value = mineigK(c-A'*y1,K);

    if value > 0
        bound.val = b'*y; bound.sol = y; % Set the current dual feasible solution
    end

    StockD.Sol = ppinfo.y;
    StockD.err4 = max( -value/(1+max(c)), 0);
end

if d_min2 > 0
    if isempty(ppinfo.z)
        ppinfo.z = c-A.'*y;
        ppinfo.y = y;
        StockD.Sol = ppinfo.y;
        StockD.err4 = 0;
    end
    
    if isempty(bound.val)
        bound.val = b'*y; bound.sol = y; % Set the current dual feasible solution
    else
        if b'*y >= bound.val, bound.val = b'*y; bound.sol = y; end
    end
end



end

