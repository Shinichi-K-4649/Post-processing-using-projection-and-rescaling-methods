function [bound,TmpSol,info] = ProposedAlgorithm_type2(A,b,c,K,sol,para)
% This function post-processes the input solution in the following steps.
% 1. Apply Algorithm 4 (Call "PostProcessingAlg" to obtain a highly accurate optimal primal solution.)
% 2. Apply Algorithm 5 (Call "Alt_PostProcessingAlg" to obtain a highly accurate optimal dual solution.)

if isfield(para,'TimeLimit')
    para.Maxtime = para.TimeLimit.P;
end

start = tic;
[bound,TmpSol,P_info] = PostProcessingAlg(A,b,c,K,sol,para);
Ptime = toc(start);
P_info.Time = Ptime;

if P_info.msg == "Primal Reducing direction" || P_info.msg == "Dual improving ray"
    info.Primal = P_info;
else
    if ~isempty(P_info.Psol.Sol), para.StockP = P_info.Psol; end
    if ~isempty(P_info.Dsol.Sol), para.StockD = P_info.Dsol; end
    if ~isempty(bound.val), para.LowerBound = bound; end

    if isfield(para,'TimeLimit')
        para.Maxtime = para.TimeLimit.D; 
    end
    
    if para.Maxtime > 0
        start = tic;
        [bound,TmpSol,D_info] = Alt_PostProcessingAlg(A,b,c,K,TmpSol,para);
        Dtime = toc(start);

        D_info.Time = Dtime;
        info.Dual = D_info;
    end
    info.Primal = P_info;
end

info.Time = Ptime+Dtime;
end

