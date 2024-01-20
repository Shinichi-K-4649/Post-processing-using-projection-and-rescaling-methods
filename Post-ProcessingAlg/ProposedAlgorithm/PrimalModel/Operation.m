function [solinfo] = Operation(A,b,c,K,InputVal,ppinfo,itr)

if ~isempty(ppinfo.x)
    % Use scaling with an approximate optimal interior feasible solution
    [solinfo] = InputPrimal(A,b,c,K,InputVal,ppinfo,itr); 
else
    % Construct the input feasibility problem
    [A_aft,K_aft] = transform_type1(A,b,c,K,InputVal);
    fprintf("\n %g : PR-Method with P-Model ... ", itr)

    % Call the projection and rescaling method
    [solinfo,err_flag] = PRwithPrimalModel(A_aft,K_aft,ppinfo);
    if err_flag == 0
        if solinfo.type == "primal"
            solinfo = ppsol2origPsol(solinfo);
        elseif solinfo.type == "dual"
            solinfo = ppsol2origDsol(solinfo,A,c,b,K);
        end
    elseif err_flag ==2 
        fprintf(" ==> eps acc stop !!!! ")
    end
end



end

