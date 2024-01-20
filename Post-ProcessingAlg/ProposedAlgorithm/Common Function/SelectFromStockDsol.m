function [y_opt] = SelectFromStockDsol(A,b,c,K,y,StockD)
% Function to extract the best dual solution among the approximate feasible dual solutions obtained using the projection and rescaling method

Stock_Dsol = StockD.Sol;
err4_table = StockD.err4;
y_opt = [];
value = max(0, -mineigK(c-A'*y,K)/(1+max(c)));

if ~isempty(Stock_Dsol)
    index = find(err4_table <= value);
    Stock = Stock_Dsol(:,index);
    ObjValue_table = b'*Stock;
    [~,index] = max(ObjValue_table);
    sol = Stock(:,index);
    if b'*sol >= b'*y, y_opt = sol; end
end

if isempty(y_opt)
    y_opt = y;
end


end

