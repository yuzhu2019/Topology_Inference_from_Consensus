clear;clc;close all;
load('data.mat');

S = cov(data);
[S_V,S_D] = eig(S);
[S_D,ind] = sort(diag(S_D)); % increasing order
S_V = S_V(:,ind);

eigvec1 = S_V(:,end);
eigvec2 = S_V(:,end-1);

data1 = [];
data2 = [];
for i = 1:502
    temp1 = abs(dot(data(i,:),eigvec1));
    temp2 = abs(dot(data(i,:),eigvec2));
    if temp1 > temp2
        data1 = [data1;data(i,:)];
    else
        data2 = [data2;data(i,:)];
    end
end
plot(sort(data1,1),'-o');
figure
plot(sort(data2,1),'-o');
