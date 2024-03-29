% compare InverseFilter, NearestCGL, StructuredGLasso
% filter: exponential decay
clear;clc;close all;

% Parameters
N = 36;                    
M_over_N = [0.75, 1, 3, 10, 30, 100, 300]; 
M = N * M_over_N; 
M_max = max(M);

beta = 0.75;

J = ones(N)/N;
H = 2*eye(N) - ones(N);

N_g = 20;
N_M = length(M);

n_reg = 15;
reg = zeros(n_reg,1);
reg1 = 0:0.01:0.11;
n_reg1 = length(reg1);

err_ipf = zeros(N_g,N_M);
err_gsi_cvx = zeros(n_reg,N_g,N_M);
err_normf = zeros(n_reg1,N_g,N_M);
%%
for n_g = 1:N_g
    % Model
    [A,~] = generate_connected_ER(N,0.1);
    [~,L] = add_weights(A,0.1,3); 
    L_f = norm(L,'fro');

    [L_V,L_D] = eig(L);
    lambda = diag(L_D); 
    lambda(lambda<=10^-10) = 0; 

    sigma = exp(-beta*lambda); % (0,1]
    h_L = L_V * diag(sigma) * L_V'; 
    h_L = 0.5*(h_L + h_L');
    
    y = h_L * randn(N,M_max); % gaussian
   
    %
    for n_M = 1:N_M
            
        M_current = M(n_M);
        S = cov(y(:,1:M_current)', 1);
        S = 0.5*(S + S');

        %% Algorithm
        [S_V,S_D] = eig(S);
        sigma_r = diag(S_D);
        sigma_r = sigma_r/max(sigma_r); % <=1  estimate the input power
        sigma_r(sigma_r<0) = 0;         % >=0  to apply the sqrt in next step
        sigma_r = sqrt(sigma_r);        % [0,1]
        lambda_r = -log(sigma_r)/beta;  % [0,Inf]
        lambda_r_inv = 1./lambda_r;
        lambda_r_inv(lambda_r==0) = 0; % pinv

        %% gsi - cvx
        S_pf = S_V * diag(lambda_r_inv) * S_V';
        S_pf = 0.5 * (S_pf + S_pf');

        s_max = max(max(abs(S_pf-diag(diag(S_pf)))));
        reg(1:14) = s_max * sqrt(log(N)/M_current) * (0.75.^(1:14));
        
        for count = 1:n_reg
            K = S_pf + reg(count) * H;
            cvx_begin quiet
               variable L_gsi(N,N) symmetric
               minimize trace(L_gsi*(K+J))-log_det(L_gsi+J)
               subject to
                   for i = 1:(N-1)
                       for j = (i+1):N
                           L_gsi(i,j) <= 0;
                       end
                   end
                   L_gsi * ones(N,1) == zeros(N,1);             
            cvx_end
            err_gsi_cvx(count,n_g,n_M) = norm(L_gsi-L,'fro')/L_f;
        end

        %% ipf
        lambda_r(lambda_r==Inf) = max(lambda_r(lambda_r<Inf));
        L_ipf = S_V * diag(lambda_r) * S_V'; 
        L_ipf = 0.5 * (L_ipf + L_ipf');
        err_ipf(n_g,n_M) = norm(L_ipf-L,'fro')/L_f;
        
        %% norm_f
        for count = 1:n_reg1
            cvx_begin quiet
              variable L_e(N,N) symmetric
              minimize norm(L_e-L_ipf,'fro') + reg1(count)*trace(L_e*H)
              subject to
                  for i = 1:(N-1)
                      for j = (i+1):N
                          L_e(i,j) <= 0;
                      end
                  end
                  L_e * ones(N,1) == zeros(N,1);             
            cvx_end
            err_normf(count,n_g,n_M) = norm(L_e-L,'fro')/L_f;
        end

    end

end



 












