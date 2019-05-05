clear;clc;close all;
load('data.mat');

S = cov(data,1);

[S_V,S_D] = eig(S);
[S_D,ind] = sort(diag(S_D)); % increasing order
S_V = S_V(:,ind);

epsil_min = 0;
epsil_max = 1;
del = 1;
max_iter_search = 5;
max_iter = 3;

%[L_all,epsil] = OSpecTemp1(S_V,epsil_min,epsil_max,del,max_iter_search,max_iter);
%L_rec = L_all(:,:,end);
load('L_rec.mat');
%% compare the leading eigenvec of S and second smallest eigenvec of recovered Laplacian
[L_V,L_D] = eig(L_rec);
[L_D,ind] = sort(diag(L_D)); % increasing order
L_V = L_V(:,ind);

S_eigvec_leading = S_V(:,end);
%L_rec_eigvec_2 = L_V(:,2);
%plot(S_eigvec_leading,'-o');
%hold on
%plot(L_rec_eigvec_2,'-*');

%%
N = size(S);
offDiag = ((ones(N)-eye(N))==1);
A = zeros(N); 
A(offDiag) = -L_rec(offDiag);
A(A<1e-07) = 0;
N_edges = length(find(A>0))/2; % 391 edges

% 3 types of states
type1 = find(S_eigvec_leading>0.1); % 16
type2 = find(S_eigvec_leading<0.1 & S_eigvec_leading>-0.1); % 14
type3 = find(S_eigvec_leading<-0.1); % 20
%type1= [1     3     5  6   7     8     9    13    18    23    34    43    46    47    48    50];
%type2= [2     4     10    11    12    14    15    19    21    26    37    39    41    42];
%type3= [16    17    20    22    24    25    27    28    29    30    31    32    33    35    36    38    40    44    45    49];

% edge_type -- 1:(1,1) 2:(1,2) 3:(1,3) 4:(2,2) 5:(2,3) 6:(3,3)

edge_list = zeros(N_edges,2); % first col - edge type; second col - edge weight
k = 1;
for i = 1:(N-1)
    for j = (i+1):N
        if A(i,j)>0
            if ismember(i,type1) && ismember(j,type1)
                edge_list(k,1) = 1; 
            elseif ismember(i,type2) && ismember(j,type2)
                edge_list(k,1) = 4; 
            elseif ismember(i,type3) && ismember(j,type3)
                edge_list(k,1) = 6;
            elseif (ismember(i,type1) && ismember(j,type2)) || (ismember(i,type2) && ismember(j,type1))
                edge_list(k,1) = 2;
            elseif (ismember(i,type1) && ismember(j,type3)) || (ismember(i,type3) && ismember(j,type1))
                edge_list(k,1) = 3;
            else
                edge_list(k,1) = 5;
            end
            edge_list(k,2) = A(i,j);
            k = k + 1;
        end
    end
end

[~,ind] = sort(edge_list(:,2),'descend');
edge_sort = edge_list(:,1);
edge_sort = edge_sort(ind);

top_k_vec = 10:20:390;
N_x = length(top_k_vec);
count = zeros(N_x,6);
for i = 1:N_x
    top_k = top_k_vec(i);
    temp = edge_sort(1:top_k);
    for j = 1:6
        count(i,j) = length(find(temp==j));
    end
end





