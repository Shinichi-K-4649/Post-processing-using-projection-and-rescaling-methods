function [bound,TmpSol,info] = Alt_PostProcessingAlg(A,b,c,K,sol,para)
% Post-Processing Algorithm using the projection and rescaling algorithm with the Dual model

%% Initialization
start = tic;
% Check the parameters
if isfield(para,'Maxtime')   Maxtime = para.Maxtime;       else Maxtime = inf;  end
if isfield(para,'Tol')       tol = para.Tol;               else tol = 1e-12;    end
if isfield(para,'PraTerCon') PraTerCon = para.PraTerCon;   else PraTerCon = 30; end

if isfield(para,'LowerBound') bound = para.LowerBound; skip =1; else skip =0; end
if skip == 0
    [ppinfo,bound,StockD] = Alt_CheckInputSol(A,b,c,K,sol);
else
   ppinfo.y = bound.sol; ppinfo.z = c-A'*bound.sol; StockD.Sol =[]; StockD.err4 = [];
end


UB = inf; LB = -inf;
StockP.Sol = []; StockP.err2 = [];
ppinfo.Q.d = [];
if ~isempty(bound.val), LB = bound.val; end
if ~isempty(ppinfo.z), tmp_val = b'*ppinfo.y; end

if isfield(para,'StockP')
    StockP = para.StockP;
else
    StockP.Sol = []; StockP.err2 = [];
end

if isfield(para,'StockD')
    Stock = para.StockD;
    StockD.Sol = [StockD.Sol,Stock.Sol];
    StockD.err4 = [StockD.err4, Stock.err4];
end

% Choose the initial input value
InputVal = (c.'*sol.x + b'*sol.y)/2;
if ~isempty(bound.val) && bound.val >= InputVal, InputVal = bound.val + abs(bound.val)*1e-8; end

sign = 0; itr = 0; unstable_flag = 0;
UB_eps = []; eps_count=0; unstable_count=0;

%% Start post-processing using the projection and rescaling algorithm
while true
    if UB - LB <= tol, break; end
    if UB == InputVal || LB == InputVal, break; end
    if toc(start) >= Maxtime
        sign =1; fprintf("\n !!! Time Up !!!! \n")
        break
    end

    % Call the projection and rescaling method
    [solinfo] = Alt_Operation(A,b,c,K,InputVal,ppinfo,itr);
    itr = itr+1;

    % Process according to the output from the projection and rescaling method
    if solinfo.type == "primal" % Output ---> Primal solution
        mineig = mineigK(solinfo.origsol,K);
        if norm(A*solinfo.origsol-b)<= 1e-4 && mineig >= -1e-4
            unstable_flag =0;
            unstable_count = 0;
            fprintf(" ==> P-Sol ::: c'x = %g ::: |Ax-b| = %g ::: mineig %g",c'*solinfo.origsol, norm(A*solinfo.origsol-b),mineig)
            eps_count = 0;
        else
            unstable_flag = 1;
            unstable_count = unstable_count+1;
            fprintf(" ==> P-Sol ::: Unstable output ")
        end

        StockP.Sol = [StockP.Sol,solinfo.origsol];
        if  mineig < 0
            err2 = -mineig/(max(b)+1);
        else
            err2 = 0;
        end
        StockP.err2 = [StockP.err2,err2];
        UB = InputVal;

    elseif solinfo.type == "dual" % Output ---> Dual solution
        mineig = mineigK(solinfo.origsol.z,K);
        if mineig >= -1e-4
            LB = InputVal;
            unstable_flag  =0;
            unstable_count =0;
            eps_count= 0;
            fprintf(" ==> D-Sol ::: b'y = %g ::: mineig %g", b'*solinfo.origsol.y, mineig)
        else
            unstable_flag  = 1;
            unstable_count = unstable_count+1;
            if LB ~= -inf
                LB = InputVal;
            end
            fprintf(" ==> D-Sol ::: Unstable output")
        end

        if mineig >= 0
            if ~isempty(bound.val)
                [bound] = UpdateLB_FeasibleSol(A,b,c,K,bound,solinfo.origsol.y);
            else
                bound.val = b.'*solinfo.origsol.y; bound.sol = solinfo.origsol.y;
                if ~isempty(StockD.Sol)
                    tmp_index = find(StockD.err4 <= 1e-4);
                    if ~isempty(tmp_index)
                        for t = 1:length(tmp_index)
                            tmp = tmp_index(t);
                            tmp_dualsol = StockD.Sol(:,tmp);
                            [bound] = UpdateLB_InfeasibleSol(A,b,c,K,bound,tmp_dualsol);
                        end
                    end
                end
            end
            err4 = 0;
        else
            if ~isempty(bound.val)
                [bound] = UpdateLB_InfeasibleSol(A,b,c,K,bound,solinfo.origsol.y);
            end
            err4 = max( -mineig/(1+max(c)), 0);
        end

        StockD.Sol = [StockD.Sol,solinfo.origsol.y];
        StockD.err4 = [StockD.err4, err4];

        if ~isempty(bound.val)
            LB = max(LB, bound.val);
        end

        if isempty(ppinfo.z)
            if err_flag == 0
                ppinfo.y = solinfo.origsol.y;
                ppinfo.z = c-A'*ppinfo.y;
                tmp_val = b'*ppinfo.y;
            end
        else
            if UB - LB <= 1
                % preserve the scaling information
                if ~isfield(ppinfo,"IntSolScale"), ppinfo.IntSolScale = solinfo.IntSolScale; end
                ppinfo.Q.d = solinfo.Q;
            else
                if mineig > 0
                    if b'*bound.sol > tmp_val
                        ppinfo.y = bound.sol;
                        ppinfo.z = c-A'*ppinfo.y;
                        tmp_val = b'*ppinfo.y;
                        if isfield(ppinfo,"IntSolScale"), ppinfo = rmfield(ppinfo,"IntSolScale"); end
                    end
                end
            end
        end
    elseif solinfo.type == "eps_acc" % Output ---> A certificate that there is no epsilon feasible solution to the input feasibility problem
        eps_count= eps_count + 1;
        if isempty(UB_eps)
            UB_eps = InputVal;
        end
        UB = InputVal;
    elseif solinfo.type == "dual reducing direction" % Output ---> A reducing direction for the dual problem
        fprintf(" ==> %s", solinfo.type)
        sign = 2;
        break
    elseif solinfo.type == "primal improving ray" % Output ---> An improving ray of the primal problem
        fprintf(" ==> %s", solinfo.type)
        sign = 5;
        break
    elseif solinfo.type == "unstable output primal" % Output ---> A primal solution affected by numerical error
        unstable_flag = 1;
        unstable_count = unstable_count+1;
        if UB ~= inf
            UB = InputVal;
        end
        fprintf(" ==> P-Sol ::: Unstable output ")
    elseif solinfo.type == "unstable output dual" % Output ---> A dual solution affected by numerical error
        unstable_flag =1;
        unstable_count = unstable_count+1;
        if LB ~= -inf
            LB = InputVal;
        end
        fprintf(" ==> D-Sol ::: Unstable output")
    end

    % Choose the input value for the next iteration (and check whether the practical termination conditions are met)
    if UB == inf
        if InputVal < LB, InputVal = LB; end
        if itr <= 20
            InputVal = InputVal + abs(InputVal)*1e-6;
        else
            InputVal = InputVal + max( 1, abs(InputVal)*1e-6);
        end

        if unstable_count >= PraTerCon
            sign = 3;
            break
        end

    elseif LB == -inf

        if unstable_flag == 0
            if itr <= 20
                InputVal = InputVal - abs(InputVal)*1e-6;
            else
                InputVal = InputVal - max( 1, abs(InputVal)*1e-6);
            end
        else
            if unstable_count <= PraTerCon
                InputVal = InputVal - max(abs(InputVal)*1e-4,1);
            else
                sign =3;
                break
            end
        end

        if eps_count >= PraTerCon
            sign = 4;
            break
        end

    else
        if solinfo.type == "eps_acc"
            InputVal = UB/4 + 3*LB/4;
        else
            if tmp_val <= LB || tmp_val >= UB
                if abs(tmp_val-UB) < UB-LB
                    InputVal = UB - abs(tmp_val-UB);
                    if InputVal == LB, InputVal = UB/2 + LB/2; end
                else
                    if unstable_count > 0
                        InputVal = UB/4 + 3*LB/4;
                    else
                        InputVal = UB/2 + LB/2;
                    end
                end
            else
                InputVal = tmp_val;
            end
        end

        if unstable_count > PraTerCon
            sign =3;
            break
        end
    end

end

%% Prepare the output objects
info.Psol = StockP;
info.Dsol = StockD;

if ~isempty(bound.sol)
    TmpSol.y = bound.sol; % Output the current dual feasible solution
else
    [y_opt] = SelectFromStockDsol(A,b,c,K,sol.y,StockD); % Extract the best dual solution among the approximate feasible dual solutions
    TmpSol.y = y_opt; % Output the approximate feasible dual solution
end

if ~isempty(StockP.Sol)
    StockP.Sol = [sol.x,StockP.Sol];
    mineig = mineigK(sol.x,K);
    if mineig < 0
        err2 = -mineig/(max(b)+1);
    else
        err2 = 0;
    end
    StockP.err2 = [err2, StockP.err2];
    info.Psol = StockP;
    [solution] = SelectFromStockPsol(A,b,c,StockP,TmpSol.y); % Extract the best primal solution among the approximate feasible primal solutions
    TmpSol.x = solution; TmpSol.z = c-A'*TmpSol.y;
else
    TmpSol.x = sol.x; TmpSol.z = c-A'*TmpSol.y;
end

if sign == 0
    info.msg = "Complete";
elseif sign == 1
    info.msg = "Time Over";
elseif sign == 2
    info.msg = "Dual Reducing direction";
    info.rd = solinfo.origsol.y;
elseif sign ==3
    info.msg = "Terminate due to numerical error";
elseif sign == 4
    info.msg = "Dual problem might have no interior feasible solution";
elseif sign == 5
    info.msg = "Primal improving ray";
    info.ir = solinfo.origsol.y;
end

info.UB = UB; info.LB = LB; info.UB_eps = UB_eps;
end

