clear;clc;%close all;

% Parameters
N = 36;                    
M_over_N = [0.75, 1, 3, 10, 30, 100, 300]; 
M = N * M_over_N; 
M_max = max(M);
p = 0.1;

beta = 0.75;

N_g = 500;
N_M = length(M);
err_ipf = zeros(N_g,N_M);
err_inv = zeros(N_g,N_M);
%%
for n_g = 1:N_g
    % Model
    [A,~] = generate_connected_ER(N,p);
    [~,L] = add_weights(A,0.1,3); 
    L_f = norm(L,'fro');

    [L_V,L_D] = eig(L);
    lambda = diag(L_D); 
    lambda(lambda<=10^-10) = 0; 
    
    % L inverse
    L_pinv = L_V * pinv(diag(lambda)) * L_V';
    L_pinv = 0.5*(L_pinv + L_pinv');
    L_pinv_f = norm(L_pinv,'fro');
   % lambda_inv = diag(pinv(diag(lambda)));
    
    sigma = exp(-beta*lambda); % (0,1], ~ 10^-5
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
        
        %
        S_pf = S_V * diag(lambda_r_inv) * S_V';
        S_pf = 0.5 * (S_pf + S_pf');
        err_inv(n_g,n_M) = norm(S_pf-L_pinv,'fro')/L_pinv_f;
        
        %%
        lambda_r(lambda_r==Inf) = max(lambda_r(lambda_r<Inf));
        L_ipf = S_V * diag(lambda_r) * S_V'; 
        L_ipf = 0.5 * (L_ipf + L_ipf');
        err_ipf(n_g,n_M) = norm(L_ipf-L,'fro')/L_f;

    end
end



 












