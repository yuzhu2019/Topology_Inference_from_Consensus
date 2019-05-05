load('L_all.mat');
L = L_all(:,:,4);

A = diag(diag(L)) - L;
w_vec = sort(A(:),'descend');
% 100 130 150 200 250 300
top_k = 300;
thr = w_vec(top_k*2);
A(A<=thr) = 0;
G = graph(A);
LWidths = 3*G.Edges.Weight/max(G.Edges.Weight);
p = plot(G,'Layout','force','LineWidth',LWidths);
color = zeros(50,1);
a=[ 1     3     5     6     7     8     9    13    18    23    34    43    46    47    48    50];
b=[ 2     4    10    11    12    14    15    19    21    26    37    39    41    42];
c=[ 16    17    20    22    24    25    27    28    29    30    31    32    33    35    36    38    40    44    45    49];
color(a) = 1;
color(b) = 3;
color(c) = 2;
p.NodeCData =color;