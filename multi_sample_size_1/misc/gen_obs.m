load('G_smallworld');
N = 36;               % graph size
obs_max = N * 1000;   % largest sample size
N_g = 20;             % number of trials

Y_mat = zeros(N,obs_max,N_g);

for n_g = 1:N_g
    L = G_mat(:,:,n_g); 
    lam_max = eigs(L,1);
    for n_obs = 1:obs_max
        T_i = randi([3 5]); % chosen from [3,4,5]
        H_i = eye(N);
        for tt = 1:T_i
            H_i = H_i * (eye(N) - rand()/lam_max*L);
        end
        Y_mat(:,n_obs,n_g) = H_i * randn(N,1);
    end
end
