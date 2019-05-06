% compute F-score metric
% L - true CGL; L_rec - recovered CGL; thr - threshold value 
% fs - return the computed F-score
% Yu, 05/03/2019
function fs = comp_fs(L, L_rec, thr)

N = size(L,1); % graph size

A = diag(diag(L)) - L;
A_rec = diag(diag(L_rec)) - L_rec;
A_rec = 0.5 * (A_rec + A_rec'); % ensure it is symmetric
A_rec = (A_rec>thr); % binary

tp = 0; % true  -  positive
fp = 0; % false -  positive
fn = 0; % false -  negative

for row = 1:N-1
    for col = (row+1):N
        if A(row,col)==1 && A_rec(row,col)==1
            tp = tp+1;
        elseif A(row,col)==0 && A_rec(row,col)==1
            fp = fp+1;
        elseif A(row,col)==1 && A_rec(row,col)==0
            fn = fn+1;
        end
    end
end

fs = 2*tp/(2*tp+fn+fp);

end


