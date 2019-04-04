%err_ipf = zeros(N_g,N_M);
%err_gsi_cvx = zeros(n_reg,N_g,N_M);
%err_normf = zeros(n_reg1,N_g,N_M);

%plot(mean(err_ipf));
hold on 
 
%%
temp1 = err_normf(1,:,:);
temp1 = reshape(temp1,20,7);
plot(mean(temp1));

temp2 = min(err_normf);
temp2 = reshape(temp2,20,7);
plot(mean(temp2));


%%
temp1 = err_gsi_cvx(15,:,:);
temp1 = reshape(temp1,20,7);
plot(mean(temp1));

temp2 = min(err_gsi_cvx);
temp2 = reshape(temp2,20,7);
plot(mean(temp2));
