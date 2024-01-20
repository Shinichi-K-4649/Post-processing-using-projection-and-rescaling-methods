function [bound,TmpSol,info] = ProposedAlgorithm(A,b,c,K,sol,para)
% This function post-processes the input solution in the following steps.
% 1. Apply Algorithm 5 (Call "Alt_PostProcessingAlg" to obtain a highly accurate optimal dual solution.)
% 2. Apply Algorithm 4 (Call "PostProcessingAlg" to obtain a highly accurate optimal primal solution.)

if isfield(para,'TimeLimit')
    para.Maxtime = para.TimeLimit.D;
end

start = tic;
[bound,TmpSol,D_info] = Alt_PostProcessingAlg(A,b,c,K,sol,para);
Dtime = toc(start);
D_info.Time = Dtime;

if D_info.msg == "Dual Reducing direction" || D_info.msg == "Primal improving ray"
    info.Dual = D_info;
else
    if ~isempty(D_info.Psol.Sol), para.StockP = D_info.Psol; end
    if ~isempty(D_info.Dsol.Sol), para.StockD = D_info.Dsol; end
    if ~isempty(bound.val), para.LowerBound = bound; end

    if isfield(para,'TimeLimit')
        para.Maxtime = para.TimeLimit.P; 
    end
    
    if para.Maxtime > 0
        start = tic;
        [bound,TmpSol,P_info] = PostProcessingAlg(A,b,c,K,TmpSol,para);
        Ptime = toc(start);

        P_info.Time = Ptime;
        info.Primal = P_info;
    end
    info.Dual = D_info;
end

info.Time = Ptime+Dtime;

end

