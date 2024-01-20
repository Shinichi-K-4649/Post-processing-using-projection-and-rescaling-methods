function [solinfo] = InputDual(A,b,c,K,InputValue,ppinfo,itr)
%% Scale the problem using an approximate optimal interior feasible solution
[m,d] = size(A); A_aft = sparse(m,d); c_aft = sparse(d,1);
if isfield(K,'l') == 1
    num_l = K.l;
    A_l = A(:,1:num_l); c_l = c(1:num_l);

    tmp = ppinfo.z(1:num_l); index = find(tmp==0, 1);
    if  ~isempty(index)
        tmp(tmp==0)=1e-16;
    end

    A_aft(:,1:num_l) = A_l./tmp.';
    c_aft(1:num_l) = c_l./tmp; 
else
    num_l = 0;
end

if ~isfield(ppinfo,"IntSolScale")
    start = num_l; start2 = num_l;
    [eig_values, eigvec_cell] = eigK_with_vec(full(ppinfo.z),K);
    for i = 1:length(K.s)
        n = K.s(i); P = eigvec_cell{i};
        A_tmp = A(:,start+1:start+n^2); c_tmp = c(start+1:start+n^2);

        d = eig_values(start2+1:start2+n);
        if min(d) < 0
            d(d<0) = 1e-16;
            eig_values(start2+1:start2+n) = d;
        end
        d = sqrt(d); d = 1./d;

        A_tmp = A_tmp*kron(P,P); c_tmp = (c_tmp.'*kron(P,P)).';
        A_tmp = A_tmp*kron(diag(d),diag(d)); c_tmp = (c_tmp.'*kron(diag(d),diag(d))).';

        A_aft(:,start+1:start+n^2) = A_tmp; c_aft(start+1:start+n^2) = c_tmp;
        start = start+n^2; start2 = start2 + n;
    end

    IntSolScale.A = A_aft(:,num_l+1:end);
    IntSolScale.c = c_aft(num_l+1:end);
    IntSolScale.EigVal = eig_values;
    IntSolScale.EigVec = eigvec_cell;
else
    A_aft(:,num_l+1:end) = ppinfo.IntSolScale.A;
    c_aft(num_l+1:end) = ppinfo.IntSolScale.c;
    eig_values = ppinfo.IntSolScale.EigVal;
    eigvec_cell = ppinfo.IntSolScale.EigVec;
end

%% Call the projection and rescaling method
[A_aft,K_aft] = transform_type1(A_aft,b,c_aft,K,InputValue); % Construct the input feasibility problem
fprintf("\n %g : PR-Method with D-Model ... ", itr)
[solinfo,err_flag] = PRwithDualModel(A_aft,K_aft,ppinfo); % Call the projection and rescaling method
if ~isfield(ppinfo,"IntSolScale"), solinfo.IntSolScale = IntSolScale; end


if err_flag == 0
    if solinfo.type == "primal"
        solinfo = Alt_ppsol2origPsol(solinfo,c,K);
    elseif solinfo.type == "dual"
        solinfo = Alt_ppsol2origDsol(solinfo,A,c);
    end

    % Scale the obtained solution back to the solution to the original problem
    if solinfo.type == "primal"
        [d,~]=size(ppinfo.z); sol = zeros(d,1);
        if isfield(K,'l') == 1
            sol(1:num_l) = solinfo.origsol(1:num_l)./ppinfo.z(1:num_l);
        else
            num_l = 0;
        end

        start = num_l; start2 = num_l;
        for i = 1:length(K.s)
            n = K.s(i); P = eigvec_cell{i};
            d = eig_values(start2+1:start2+n);
            d = sqrt(d); d = 1./d;
            sol_tmp = solinfo.origsol(start+1:start+n^2);
            sol_tmp = kron(diag(d),diag(d))*sol_tmp;
            sol_tmp = kron(P,P)*sol_tmp;
            sol(start+1:start+n^2) = sol_tmp;
            start = start+n^2; start2 = start2 +n;
        end

        solinfo.origsol = sol;
    end

elseif err_flag ==2
    fprintf(" ==> eps acc stop !!!! ")
end

end

