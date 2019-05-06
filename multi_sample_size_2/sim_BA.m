clear;clc;
addpath('misc');
addpath('reweight')
load('G_BA.mat'); % BA 36x36x20
load('Y_BA.mat'); % 36x36000x20

N = 36; 
obs_vec = N * [1000 300 100 30 10 3 1];
obs_len = length(obs_vec);
trials = 20;
max_iter = 1; % iteration number of reweighted scheme
del = 1;
thr = 0.5;

err = zeros(max_iter+1,obs_len,trials);
fs = zeros(max_iter+1,obs_len,trials);
epsi_1_list = zeros(obs_len,trials);

for n_trials = 1:trials
    sprintf('n_trials = %d',n_trials)
    L = G_mat(:,:,n_trials);
    sum_eig = trace(L);
    L_f = norm(L,'fro');
    Y = Y_mat(:,:,n_trials);
    epsi_1_init = 0.01; 
    for n_obs = 1:obs_len
        sprintf('n_obs = %d',n_obs)
        obs_cur = obs_vec(n_obs);
        S = Y(:,1:obs_cur)*Y(:,1:obs_cur)'/obs_cur;
        [S_V,S_D] = eig(S);
        [~,ind] = sort(diag(S_D)); 
        V_tilde = S_V(:,ind);
        [L_rec,epsi_1_rt] = OSpecTemp1(V_tilde, epsi_1_init, max_iter, del);
        epsi_1_init = epsi_1_rt;
        epsi_1_list(n_obs,n_trials) = epsi_1_rt;
        for tt = 1:(max_iter+1)
            temp = L_rec(:,:,tt);
            temp = temp*sum_eig/trace(temp);
            err(tt,n_obs,n_trials) = norm(temp-L,'fro')/L_f;
            fs(tt,n_obs,n_trials) = comp_fs(L, temp, thr);
        end
    end
end


      




