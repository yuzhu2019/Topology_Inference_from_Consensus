load('ER_err.mat'); % 2x7x20
err1 = err;
load('SBM_err.mat');
err2 = err;
load('smallworld_err.mat');
err3 = err;
load('BA_err.mat');
err4 = err;

err1_avg = flip(mean(squeeze(err1(1,:,:)),2));
err2_avg = flip(mean(squeeze(err2(1,:,:)),2));
err3_avg = flip(mean(squeeze(err3(1,:,:)),2));
err4_avg = flip(mean(squeeze(err4(1,:,:)),2));

plot(err1_avg,'-o');
hold on
plot(err2_avg,'-o');
plot(err3_avg,'-o');
plot(err4_avg,'-o');


