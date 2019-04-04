clear;clc;%close all;

% Parameters
N = 36;                    
M_over_N = [0.75, 1, 3, 10, 30, 100, 300]; 
M = N * M_over_N; 
M_max = max(M);
p = 0.1;
T = 3;
alpha_ratio = [0.7,0.8,0.9];

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
    lambda_max = max(lambda);
    
    %
    L_pinv = L_V * pinv(diag(lambda)) * L_V';
    L_pinv = 0.5*(L_pinv + L_pinv');
    L_pinv_f = norm(L_pinv,'fro');
 %   lambda_inv = diag(pinv(diag(lambda)));
    
    alpha = alpha_ratio/lambda_max;
    sigma = (1-alpha(1)*lambda).*(1-alpha(2)*lambda).*(1-alpha(3)*lambda); %0<sigma<=1
    h_L = L_V * diag(sigma) * L_V'; 
    h_L = 0.5*(h_L + h_L');
    
    y = h_L * randn(N,M_max); % gaussian
    
    %
	f_coeff = [-alpha',ones(T,1)];
	f = 1;
	for t = 1:T
        f = conv(f,f_coeff(t,:));
    end
    
	lambda_r = zeros(N,1);
	alpha_max_inv = 1/max(alpha);

    for n_M = 1:N_M

        M_current = M(n_M);
        S = cov(y(:,1:M_current)', 1);
        S = 0.5*(S + S');

        %% Algorithm
        [S_V,S_D] = eig(S);
        sigma2_r = diag(S_D);
        sigma2_r = sigma2_r/max(sigma2_r);
        sigma2_r(sigma2_r<0) = 0;
        sigma_r = sigma2_r.^(0.5);

        for n = 1:N
            temp = f;
            temp(end) = temp(end)-sigma_r(n);
            temp = roots(temp);
            temp2 = temp(imag(temp)==0 & real(temp)<=(alpha_max_inv+10^-5) );
            assert(length(temp2)==1,'Error!');
            lambda_r(n) = temp2;
        end

        % ipf
        L_ipf = S_V * diag(lambda_r) * S_V';       
        err_ipf(n_g,n_M) = norm(L_ipf-L,'fro')/L_f;

        %% gsi - cvx
        S_pf = S_V * pinv(diag(lambda_r)) * S_V';
        S_pf = 0.5 * (S_pf + S_pf');
        err_inv(n_g,n_M) = norm(S_pf-L_pinv,'fro')/L_pinv_f;
             
    end
end



 












