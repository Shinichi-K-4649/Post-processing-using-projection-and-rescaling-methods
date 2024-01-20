function [judge,solution,upper,Rescale,Q_orig_d,Q_orig_p,A] = PP_BP_SP(K,Projection,upper,eps,Rescale,Q_orig_d,Q_orig_p,A_orig,A)
% Basic procedure using the smooth perceptron scheme

solution = []; judge = 0;
mu = 2;
u_bar = eyeK(K)/(K.l+sum(K.s));
u = u_bar;
k=0;
y = smooth_operator(u-Projection*u,mu,u_bar,K);
v = Projection*y; z = y-v;

while true
    minimum_eig_value = mineigK(z,K);

    % Check whether z is a feasible solution to the input problem
    if minimum_eig_value > -1e-16
        judge = 1; solution = z; break
    end

    % Check whether v is an infeasibility certificate to the input problem
    minimum_eig_value = mineigK(v,K);
    if minimum_eig_value >= -1e-16 
        judge = -1; solution = v; break
    end

    % Check whether cuts can be obtained for the input problem
    value = v.'*eyeK(K);
    [eig_values, eigvec_cell] = eigK_with_vec(v,K);
    if value >0
        index = find(eig_values>0);
    else
        index = find(eig_values<0);
    end

    cut_value = ones(length(eig_values),1);
    for i = 1:length(index)
        num = index(i);
        check = eig_values/(-eig_values(num)); check(check<0)=0;
        check_value = sum(check);
        cut_value(num) = check_value;
    end

    Index = find(cut_value<=eps);
    if isempty(Index)
        if k>=100
            Index = find(cut_value < 1);
        end
    end

    if ~isempty(Index)
        % Scaling the coefficient matrix using the cuts
        ans_vec = ones(K.l,1); 
        index = find(Index <= K.l); 
        if ~isempty(index)
            ans_vec(Index(index)) = cut_value(Index(index));
            tmp = Rescale{1};
            Rescale{1} = tmp.*ans_vec;
            A(:,1:K.l) = A_orig(:,1:K.l).*Rescale{1}';

            tmp = Rescale{1};
            tmp = ones(length(K.l),1)./tmp;
            ans_vec = zeros(K.l,1); ans_vec(Index(index))=1;
            upper(1:K.l,1) = tmp.*ans_vec + upper(1:K.l,1);
        end

        num_l = K.l; start = K.l; start_s = K.l;
        for i = 1:length(K.s) 
            index1 = find(Index <= start+K.s(i));
            index2 = find(Index <= start);
            index = setdiff(index1,index2);
            if ~isempty(index)
                Q = eigvec_cell{i};
                ans_vec = ones(K.s(i),1);
                ans_vec(Index(index)-start) = sqrt(eps);
                Q_output = Q*diag(ans_vec)*Q.';

                ans_vec = zeros(K.s(i),1); ans_vec(Index(index)-start) = 1;
                upper(num_l+i,1) = upper(num_l+i,1)+...
                    sum(diag(Q_orig_p.left{i+1}*Q*diag(ans_vec)*Q.'*Q_orig_p.right{i+1}.'));

                ans_vec = ones(K.s(i),1);
                ans_vec(Index(index)-start)= 1/sqrt(eps);
                rescale = Q*diag(ans_vec)*Q.';

                Q_orig_d.left{i+1} = Q_orig_d.left{i+1}*Q_output.';
                Q_orig_d.right{i+1} = Q_orig_d.right{i+1}*Q_output;
                Tmp = kron(Q_orig_d.left{i+1}, Q_orig_d.right{i+1});

                Rescale{1+i} = Tmp;
                A(:,start_s+1:start_s+K.s(i)^2) = A_orig(:,start_s+1:start_s+K.s(i)^2)*Tmp;
                Q_orig_p.left{i+1} = Q_orig_p.left{i+1}*rescale;
                Q_orig_p.right{i+1} = Q_orig_p.right{i+1}*rescale.';
            end
            start_s = start_s + K.s(i)^2;
            start = start + K.s(i);
        end
        break
    end

    % Update the current point using the smooth perceptron scheme
    theta = 2/(k+3);
    u = (1-theta)*(u+theta*y) + theta^2*smooth_operator(u-Projection*u,mu,u_bar,K);
    mu = (1-theta)*mu;
    y = (1-theta)*y + theta*smooth_operator(u-Projection*u,mu,u_bar,K);
    v = Projection*y; z = y-v;
    k = k+1;
end


end

