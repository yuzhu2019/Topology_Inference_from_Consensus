clear;clc;close all;
addpath('usc_codes');
load('us_temp_data.mat'); 

S = cov(data_2d_state,1); % sample size: 5844, 45 states

[S_V,S_D] = eig(S);
[S_D,ind] = sort(diag(S_D)); % increasing order
S_V = S_V(:,ind);

epsil_min = 0;
epsil_max = 5;
epsi_2 = 0.1;
del = 1;
max_iter_search = 5;
max_iter = 5;

[L_all,epsil] = OSpecTemp2(S_V,epsil_min,epsil_max,epsi_2,del,max_iter_search,max_iter);

draw_us_temp_graph(L_all(:,:,end), center_vector);

% 334 edges (1e-03, 1e-08)
