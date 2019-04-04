plot(mean(err_ipf));
hold on 
%%
temp1 = err_normf(1,:,:);
temp1 = reshape(temp1,20,7);
plot(mean(temp1));

temp2 = [
   min(err_normf(:,:,1)); 
   min(err_normf(:,:,2)); 
   min(err_normf(:,:,3)); 
   min(err_normf(1:9,:,4)); 
   min(err_normf(1:9,:,5)); 
   min(err_normf(1:9,:,6)); 
   min(err_normf(1:9,:,7)); 
];
plot(mean(temp2,2));

%%
temp1 = err_norm2(1,:,:);
temp1 = reshape(temp1,20,7);
plot(mean(temp1));

temp2 = [
   min(err_norm2(:,:,1)); 
   min(err_norm2(:,:,2)); 
   min(err_norm2(:,:,3)); 
   min(err_norm2(1:9,:,4)); 
   min(err_norm2(1:9,:,5)); 
   min(err_norm2(1:9,:,6)); 
   min(err_norm2(1:9,:,7)); 
];
plot(mean(temp2,2));

%%
temp1 = err_gsi_cvx(15,:,:);
temp1 = reshape(temp1,20,7);
plot(mean(temp1));

temp2 = min(err_gsi_cvx);
temp2 = reshape(temp2,20,7);
plot(mean(temp2));



  