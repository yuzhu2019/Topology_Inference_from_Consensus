%err_ipf = zeros(N_g,N_M,N_a);
%err_normf = zeros(14,N_g,N_M,N_a);
%clear;clc;close all;
%load('err_ipf.mat');
%load('err_normf.mat');

plot(mean(err_ipf(:,:,1)));
hold on
plot(mean(err_ipf(:,:,2)));
plot(mean(err_ipf(:,:,3)));
plot(mean(err_ipf(:,:,4)));

for t = 1:4
    temp = err_normf(1,:,:,t);
    temp = reshape(temp,20,7);
    plot(mean(temp));
end

for t = 1:4
    temp = err_normf(:,:,:,t);
    temp2 = [ min(temp(:,:,1)); 
              min(temp(:,:,2));
              min(temp(:,:,3));
              min(temp(1:9,:,4));
              min(temp(1:9,:,5));
              min(temp(1:9,:,6));
              min(temp(1:9,:,7))]';
    plot(mean(temp2));
end


         