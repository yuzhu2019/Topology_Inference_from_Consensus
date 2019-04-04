% Generate a N*N grid graph
% Output: adjacency matrix A, combinatorial Laplacian L
% Yu Zhu, Rice ECE, 12/04/2018
function [A,L] = generate_grid(N)
    N_v = N^2; 
    A = zeros(N_v);
    for i = 1:N
        for j = 1:(N-1)
            s = (i-1)*N+j;
            A(s,s+1) = 1;
            t = (j-1)*N+i;
            A(t,t+N) = 1;
        end
    end
    A = A + A';
    L = diag(sum(A))-A;
end

