function [L_all,epsil] = OSpecTemp1(V_tilde,epsil_min,epsil_max,del,max_iter_search,max_iter)

    N = size(V_tilde, 1);
    diag_idx = 1:N+1:N^2;
    non_diag_idx = setdiff(1:N^2,diag_idx)';
    tau = 1;
    delt = 0.001;
    L_all = zeros(N,N,max_iter+1);

    %% first stage
    epsil_low = epsil_min;
    epsil_high = epsil_max;
    for tt = 1:max_iter_search
        epsil = (epsil_low+epsil_high)/2;
        cvx_begin 
            variables L(N,N)
            variables L_tilde(N,N)
            variables lambda(N,1)
            minimize norm(L(non_diag_idx),1)
            subject to
                L == semidefinite (N);
                L(non_diag_idx) <= 0;
                L * ones(N,1) == 0;
                L_tilde == V_tilde * diag(lambda) * V_tilde';
                norm(L - L_tilde) <= epsil;
                lambda(1) == 1;
                for i = 1:N-del
                    lambda(i) >= lambda(i+del);
                end
        cvx_end
        
        if sum(sum(isnan(L))) == 0 
            epsil_high = epsil;
            L_last_valid = L;
            epsil_last_valid = epsil;
        else
            epsil_low = epsil;
        end
    end
    
    L_prev = L_last_valid;
    epsil = epsil_last_valid;
    L_all(:,:,1) = L_last_valid;
    
    %% second stage
    for tt = 1:max_iter
        weigh = tau*ones(N,N)./(abs(L_prev)+delt*ones(N,N));
        cvx_begin 
            variables L(N,N)
            variables L_tilde(N,N)
            variables lambda(N,1)
            minimize (-1*weigh(non_diag_idx)'*L(non_diag_idx))
            subject to
                L == semidefinite (N);
                L(non_diag_idx) <= 0;
                L * ones(N,1) == 0;
                L_tilde == V_tilde * diag(lambda) * V_tilde';
                norm(L - L_tilde) <= epsil;
                lambda(1) == 1;
                for i = 1:N-del
                    lambda(i) >= lambda(i+del);
                end
        cvx_end
        L_prev = L;
        L_all(:,:,tt+1) = L;
    end
    
end