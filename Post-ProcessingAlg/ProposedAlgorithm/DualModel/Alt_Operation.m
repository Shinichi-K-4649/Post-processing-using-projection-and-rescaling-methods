function [solinfo] = Alt_Operation(A,b,c,K,InputVal,ppinfo,itr)

if ~isempty(ppinfo.z)
    % Use scaling with an approximate optimal interior feasible solution
    [solinfo] = InputDual(A,b,c,K,InputVal,ppinfo,itr);
else
    % Construct the input feasibility problem
    [A_aft,K_aft] = transform_type1(A,b,c,K,InputVal);
    fprintf("\n %g : PR-Method with D-Model ... ", itr)

    % Call the projection and rescaling method
    [solinfo,err_flag] = PRwithDualModel(A_aft,K_aft,ppinfo);
    if err_flag == 0
        if solinfo.type == "primal"
            solinfo = Alt_ppsol2origPsol(solinfo,c,K);
        elseif solinfo.type == "dual"
            solinfo = Alt_ppsol2origDsol(solinfo,A,c);
        end
    elseif err_flag ==2 
        fprintf(" ==> eps acc stop !!!! ")
    end
end

end

