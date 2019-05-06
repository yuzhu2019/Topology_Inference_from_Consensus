load('ER_fs.mat');
fs1 = fs;
load('SBM_fs.mat');
fs2 = fs;
load('smallworld_fs.mat');
fs3 = fs;
load('BA_fs.mat');
fs4 = fs;

fs1_avg = flip(mean(squeeze(fs1(1,:,:)),2));
fs2_avg = flip(mean(squeeze(fs2(1,:,:)),2));
fs3_avg = flip(mean(squeeze(fs3(1,:,:)),2));
fs4_avg = flip(mean(squeeze(fs4(1,:,:)),2));

plot(fs1_avg,'-o');
hold on
plot(fs2_avg,'-o');
plot(fs3_avg,'-o');
plot(fs4_avg,'-o');



