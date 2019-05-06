N = 36;
p = 0.5;
N_g = 20;
G_mat = zeros(N,N,N_g);

for n_g = 1:N_g
    [~, G_mat(:,:,n_g)] = generate_connected_ER(N,p);
end


% N = 36;
% N_g = 20;
% G_mat = zeros(N,N,N_g);
% 
% for n_g = 1:N_g
%     [~, G_mat(:,:,n_g)] = generate_connected_SBM(N,2,0.05,0.2);
% end

% N = 36;
% N_g = 20;
% G_mat = zeros(N,N,N_g);
% for n_g = 1:N_g
%     [~, G_mat(:,:,n_g)] = generate_connected_smallworld(N,2,0.2);
% end

% 
% N = 36;
% N_g = 20;
% G_mat = zeros(N,N,N_g);
% for n_g = 1:N_g
%     [seed,~] = generate_connected_ER(5,0.2);
%     [~, G_mat(:,:,n_g)] = generate_connected_BA(N, 3, seed);
% end



