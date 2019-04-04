% generate the input as time series
% compare the cases when the input has different correlations
clear;clc;close all;

% Parameters
N = 36;  p = 0.1;                   
M_over_N = [0.75, 1, 3, 10, 30, 100, 300]; 
M = N * M_over_N;  M_max = max(M);
T = 3;  alpha_ratio = [0.7,0.8,0.9];
a = [0,0.25,0.5,0.75]; % correlation parameter

H = 2*eye(N) - ones(N);

N_g = 20;
N_M = length(M);
N_a = length(a);
err_ipf = zeros(N_g,N_M,N_a);
err_normf = zeros(14,N_g,N_M,N_a);
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

    alpha = alpha_ratio/lambda_max;
    sigma = (1-alpha(1)*lambda).*(1-alpha(2)*lambda).*(1-alpha(3)*lambda); %0<sigma<=1
    h_L = L_V * diag(sigma) * L_V'; 
    h_L = 0.5*(h_L + h_L');

    alpha_max_inv = 1/max(alpha);
   
	f_coeff = [-alpha',ones(T,1)];
	f = 1;
	for t = 1:T
        f = conv(f,f_coeff(t,:));
    end
    
    %
    for n_a = 1:N_a
        
        % generate data
        x = zeros(N,M_max);
        x(:,1) = randn(N,1);
        for m = 2:M_max
            x(:,m) = a(n_a)*x(:,m-1) + (1-a(n_a))*randn(N,1);
        end
        y = h_L*x;
        
        %
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

            lambda_r = zeros(N,1);
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
            err_ipf(n_g,n_M,n_a) = norm(L_ipf-L,'fro')/L_f;

            %% norm_f
            if n_M<=3
                reg1 = [0, 0.055:0.0025:0.085]; %14
            else
                reg1 = 0:0.01:0.08; %9
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
                err_normf(count,n_g,n_M,n_a) = norm(L_e-L,'fro')/L_f;
            end

        end
    end
end



 












