function [sol] = SelectFromStockPsol(A,b,c,Stock,dsol)
% Function to extract the best primal solution among the approximate feasible primal solutions obtained using the projection and rescaling method

Stock_Psol = Stock.Sol;
err2_table = Stock.err2;

if ~isempty(Stock_Psol)
    [~,num] = size(Stock_Psol);
    y = dsol;
    z = c-A'*y;

    V = A*Stock_Psol-ones(length(b),num).*b;
    NormValue1 = sqrt(sum(V.*V));
    ObjValue1 = c.'*Stock_Psol;
    V = Stock_Psol.*(ones(length(z),num).*z);
    SlackValue1 = sum(V);
    
    err_table = NormValue1/(1+max(abs(b)));
    err_table = err_table + err2_table;
    err_table = err_table + abs((ObjValue1 - b.'*y))./(1+abs(ObjValue1)+abs(b.'*y));
    err_table = err_table + abs(SlackValue1)./(1+abs(ObjValue1)+abs(b.'*y));
    
    [~,index] = min(err_table);
    sol = Stock_Psol(:,index);
end

end

