%load('G_ER_36_01.mat');       % 3.5750
%load('G_SBM_36_2_005_02');    % 4.5222
%load('G_smallworld');         % 4
load('G_BA');                  % 5.3889

N = 36;
N_g = 20;
deg_vec = zeros(N_g,1);
for n_g = 1:N_g
    deg_vec(n_g) = trace(G_mat(:,:,n_g));
end
deg_vec = deg_vec/N;
deg_avg = mean(deg_vec);



