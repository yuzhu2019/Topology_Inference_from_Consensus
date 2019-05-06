% modified on 05/05/2019
function [L_rec,epsi_1_rt] = SpecTemp(V_tilde, epsi_1_init, max_iter)

    N = size(V_tilde, 1);
    diag_idx = 1:N+1:N^2;
    non_diag_idx = setdiff(1:N^2,diag_idx)';
    
    epsi_1_hop = 0.005; %%% star!
    tau = 1;
    delt = 0.001;
    
    L_rec = zeros(N,N,max_iter+1);
    %% stage 1
    flag = 0;
    epsi_1 = epsi_1_init;
    while flag == 0
        cvx_begin quiet
            variables L(N,N)
            variables L_tilde(N,N)
            variables lambda(N,1)
            minimize norm(L(non_diag_idx),1)
            subject to
                L == semidefinite (N);
                L(non_diag_idx) <= 0;
                L * ones(N,1) == 0;
                L_tilde == V_tilde * diag(lambda) * V_tilde';
                norm(L - L_tilde) <= epsi_1;
                lambda(1) == 1;
        cvx_end
        epsi_1 = epsi_1 + epsi_1_hop;
        flag = sum(sum(isnan(L))) == 0;
    end
    
    L_rec(:,:,1) = L;
    L_prev = L;
    epsi_1 = epsi_1 - epsi_1_hop;
    epsi_1_rt = epsi_1;
    
    %% stage 2
    for tt = 1:max_iter
        weigh = tau*ones(N,N)./(abs(L_prev)+delt*ones(N,N));
        epsi_1 = epsi_1_rt;
        flag = 0;
        while flag == 0
            cvx_begin quiet
                variables L(N,N)
                variables L_tilde(N,N)
                variables lambda(N,1)
                minimize (-1*weigh(non_diag_idx)'*L(non_diag_idx))
                subject to
                    L == semidefinite (N);
                    L(non_diag_idx) <= 0;
                    L * ones(N,1) == 0;
                    L_tilde == V_tilde * diag(lambda) * V_tilde';
                    norm(L - L_tilde) <= epsi_1;
                    lambda(1) == 1;
            cvx_end
            epsi_1 = epsi_1 + epsi_1_hop;
            flag = sum(sum(isnan(L))) == 0;
        end
        L_rec(:,:,tt+1) = L;
        L_prev = L;
    end
    
end
