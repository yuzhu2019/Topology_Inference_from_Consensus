% comparison - with regularization
% compare different graph types: ER, grid, SBM, scale-free, small world
% compare input: gaussian, uniform
clear;clc;close all;

% Parameters
N = 36;                    
M_over_N = [0.75, 1, 3, 10, 30, 100, 300]; 
M = N * M_over_N; 
M_max = max(M);
T = 3;
alpha_ratio = [0.7,0.8,0.9];

J = ones(N)/N;
H = 2*eye(N) - ones(N);

N_g = 20;
N_M = length(M);

n_reg = 15;
reg = zeros(n_reg,1);

err_ipf = zeros(N_g,N_M);
err_gsi_cvx = zeros(n_reg,N_g,N_M);
err_normf = zeros(14,N_g,N_M);
err_norm2 = zeros(14,N_g,N_M);
%%
    [A,~] = generate_grid(sqrt(N));
for n_g = 1:N_g
    % Model
    %[A,~] = generate_connected_smallworld(N,2,0.2);
    %seed = generate_connected_ER(5,0.2);
    %[A,~] = generate_connected_BA(N, 3, seed);
    %[A,~] = generate_connected_SBM(N,2,0.05,0.2);
    %[A,~] = generate_connected_ER(N,0.1);
    [~,L] = add_weights(A,0.1,3); 
    L_f = norm(L,'fro');

    [L_V,L_D] = eig(L);
    lambda = diag(L_D); 
    lambda(lambda<=10^-10) = 0; 
    lambda_max = max(lambda);

    alpha = alpha_ratio/lambda_max;
    sigma = (1-alpha(1)*lambda).*(1-alpha(2)*lambda).*(1-alpha(3)*lambda); %0<sigma<=1
    h_L = L_V * diag(sigma) * L_V'; 
    h_L = 0.5*(h_L + h_L');
    
    y = h_L * randn(N,M_max); % gaussian
    %y = h_L * sqrt(12) * (rand(N,M_max)-0.5); % uniform
    
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
        lambda_r_inv = 1./lambda_r;
        lambda_r_inv(lambda_r==0) = 0;
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

        %% norm_f
        if n_M<=3
            reg1 = [0, 0.055:0.0025:0.085];
        else
            reg1 = 0:0.01:0.08;
        end
        n_reg1 = length(reg1);
        
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

        %% norm_2
        if n_M<=3
            reg2 = [0,0.008:0.001:0.02];
        else
            reg2 =  0:0.002:0.016;
        end
        n_reg2 = length(reg2);
        
        for count = 1:n_reg2 
            cvx_begin quiet
               variable L_e2(N,N) symmetric
               minimize norm(L_e2-L_ipf) + reg2(count)*trace(L_e2*H)
               subject to
                   for i = 1:(N-1)
                       for j = (i+1):N
                           L_e2(i,j) <= 0;
                       end
                   end
                   L_e2 * ones(N,1) == zeros(N,1);             
            cvx_end
            err_norm2(count,n_g,n_M) = norm(L_e2-L,'fro')/L_f;
        end

    end

end



 












