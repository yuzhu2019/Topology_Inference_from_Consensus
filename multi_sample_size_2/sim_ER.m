clear;clc;
addpath('misc');
addpath('reweight_modified');
load('G_ER_36_01.mat'); 
load('Y_ER_36_01.mat'); 

N = 36; 
obs_vec = N * [1000 300 100 30 10 3 1];
obs_len = length(obs_vec);
trials = 20;
max_iter = 1; % iteration number of reweighted scheme
del = 1;
thr = 0.5;

err1 = zeros(max_iter+1,obs_len,trials);
fs1 = zeros(max_iter+1,obs_len,trials);
epsi_1_list1 = zeros(obs_len,trials);

err2 = zeros(max_iter+1,obs_len,trials);
fs2 = zeros(max_iter+1,obs_len,trials);
epsi_1_list2 = zeros(obs_len,trials);

for n_trials = 1:trials
    sprintf('n_trials = %d',n_trials)
    L = G_mat(:,:,n_trials);
    sum_eig = trace(L);
    L_f = norm(L,'fro');
    Y = Y_mat(:,:,n_trials);
    epsi_1_init = 0.005; 
    for n_obs = 1:obs_len
        sprintf('n_obs = %d',n_obs)
        obs_cur = obs_vec(n_obs);
        S = Y(:,1:obs_cur)*Y(:,1:obs_cur)'/obs_cur;
        [S_V,S_D] = eig(S);
        [~,ind] = sort(diag(S_D)); 
        V_tilde = S_V(:,ind);
        %
        [L_rec1,epsi_1_st] = SpecTemp(V_tilde, epsi_1_init, max_iter);
        [L_rec2,epsi_1_ost] = OSpecTemp1(V_tilde, epsi_1_st, max_iter, del);
        epsi_1_init = epsi_1_st;
        epsi_1_list1(n_obs,n_trials) = epsi_1_st;
        epsi_1_list2(n_obs,n_trials) = epsi_1_ost;
        for tt = 1:(max_iter+1)
            temp1 = L_rec1(:,:,tt);
            temp1 = temp1*sum_eig/trace(temp1);
            err1(tt,n_obs,n_trials) = norm(temp1-L,'fro')/L_f;
            fs1(tt,n_obs,n_trials) = comp_fs(L, temp1, thr);
            temp2 = L_rec2(:,:,tt);
            temp2 = temp2*sum_eig/trace(temp2);
            err2(tt,n_obs,n_trials) = norm(temp2-L,'fro')/L_f;
            fs2(tt,n_obs,n_trials) = comp_fs(L, temp2, thr);
        end
    end
end


      




