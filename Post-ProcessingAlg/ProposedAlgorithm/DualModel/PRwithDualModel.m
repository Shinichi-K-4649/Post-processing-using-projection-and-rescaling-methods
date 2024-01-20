function [solinfo,err_flag] = PRwithDualModel(A_orig,K,ppinfo)
% Projection and rescaling method with the dual model

%% Initialization
err_flag = 0; [m,~] = size(A_orig);
u = ones(K.l+length(K.s),1);
judge_vec = [];
xi = 1/4;
epsilon = 1e-16;
scal = 0;
judge_vec(1:K.l,1)= xi*(1-epsilon)/((1-xi)*epsilon);
Q=[];

if isfield(ppinfo.Q.d, "Rescale")
    % Use the scaling information from previous iterations
    Rescale = ppinfo.Q.d.Rescale;
    A = A_orig;
    A(:,1:K.l) = A_orig(:,1:K.l).*Rescale{1}';
    start = K.l;
    for i = 1:length(K.s)
        nn = K.s(i)^2;
        A(:,start+1:start+nn) = A_orig(:,start+1:start+nn)*Rescale{1+i};
        start = start + nn;
    end
    
    u = ppinfo.Q.d.u;
    rescale_check = 1;
    Q_orig_p = ppinfo.Q.d.orig.p;
    Q_orig_d = ppinfo.Q.d.orig.d;
    
    for i = 1:length(K.s)
        judge_vec(K.l+i,1)= K.s(i)*xi*(1-epsilon)/((1-xi)*epsilon);
    end

else
    A = A_orig;
    Rescale = cell(1 + length(K.s),1);
    Rescale{1} = ones(K.l,1);
    rescale_check = 0;
    Q_orig = cell(1 + length(K.s),1);
    Q_orig{1} = ones(K.l,1);
    for i = 1:length(K.s)
        Rescale{i+1} = eye(K.s(i)^2);
        Q_orig{i+1} = eye(K.s(i));
        judge_vec(K.l+i,1)= K.s(i)*xi*(1-epsilon)/((1-xi)*epsilon);
    end

    Q_orig_p.left = Q_orig;
    Q_orig_d.left = Q_orig;
    Q_orig_p.right = Q_orig;
    Q_orig_d.right = Q_orig;
end

%% Start the projection and rescaling method
while true
    % Compute a projection matrix
    [U,Sigma,V] = svd(full(A));
    if ~isreal(V)
        fprintf("\n numerical error of SVD @ P-ChuabnovAlg")
        err_flag =1; break
    end

    tol = max(max(Sigma))/1e16;
    s = diag(Sigma(1:m,1:m));
    index = find(s<tol);
    de_rank = length(index);
    Projection = V(:,1:m-de_rank)*V(:,1:m-de_rank)';

    % Call the basic procedure
    [judge,solution,u,Rescale,Q_orig_d,Q_orig_p,A] ...
        = Alt_PP_BP_SP(K,Projection,u,xi,Rescale,Q_orig_d,Q_orig_p,A_orig,A);

    if judge == 1 % Find a feasible solution to the input feasibility problem
        sign = 1; break
    elseif judge == -1 % Find a feasible solution to the alternative problem of the input feasibility problem
        sign = -1; break
    elseif judge == 0
        index = find((u-judge_vec)>=0, 1);

        if ~isempty(index)
            sign = 0; % Determine that there is no epsilon feasible solution to the input feasibility problem
            break
        end
        scal = scal+1;
    end
end

%% Prepare the output objects
fprintf(" %g scaling ... ", scal)

if err_flag == 0
    check =1;
    if scal > 0, rescale_check = 1; end
    if sign == -1
        solinfo.type = "primal"; 
        if rescale_check == 1
            tmp_sol = solution;
            tmp_sol(1:K.l) = Rescale{1}.*solution(1:K.l);
            start = K.l;
            for i = 1:length(K.s)
                nn = K.s(i)^2;
                tmp_sol(start+1:start+nn) = Rescale{i+1}*solution(start+1:start+nn);
                start = start + nn;
            end
            solinfo.sol = tmp_sol;
        else
            solinfo.sol = solution;
        end
        check = 0; 
    elseif sign == 1
        solinfo.type = "dual";
        yyy = solution;
        s = diag(Sigma(1:m-de_rank,1:m-de_rank)); s = 1./s;
        S = zeros(m,m-de_rank); S(1:m-de_rank,1:m-de_rank)=diag(s);
        y1 = U*S*V(:,1:m-de_rank).'*yyy;
        kappa1 = y1(length(y1));
        if kappa1 < 0
            solinfo.type = "numerical error output";
        else
            solinfo.sol = y1;
        end
    elseif sign == 0
        check = 0;
        solinfo.type ="eps_acc";
        err_flag = 2;
    end

    if check == 1
        if scal == 0
            solinfo.Q = ppinfo.Q.d;
        else
            % preserve the scaling information
            Q.orig.p = Q_orig_p;
            Q.orig.d = Q_orig_d;
            Q.Rescale = Rescale;
            Q.u = u;
            solinfo.Q = Q;
        end
    end

else
    solinfo.type = "SVD numerical error";
end

end

