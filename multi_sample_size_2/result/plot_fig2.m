clear;clc;
load('fs1_36_01.mat');
fs1_1 = fs1;
load('fs1_36_03.mat');
fs1_2 = fs1;
load('fs1_36_05.mat');
fs1_3 = fs1;
load('fs2_36_01.mat');
fs2_1 = fs2;
load('fs2_36_03.mat');
fs2_2 = fs2;
load('fs2_36_05.mat');
fs2_3 = fs2;

a1 = mean(squeeze(fs1_1(1,:,:)),2); a1=flip(a1);
a2 = mean(squeeze(fs1_1(2,:,:)),2); a2=flip(a2);

b1 = mean(squeeze(fs1_2(1,:,:)),2); b1=flip(b1);
b2 = mean(squeeze(fs1_2(2,:,:)),2); b2=flip(b2);

c1 = mean(squeeze(fs1_3(1,:,:)),2); c1=flip(c1);
c2 = mean(squeeze(fs1_3(2,:,:)),2); c2=flip(c2);

d1 = mean(squeeze(fs2_1(1,:,:)),2); d1=flip(d1);
d2 = mean(squeeze(fs2_1(2,:,:)),2); d2=flip(d2);

e1 = mean(squeeze(fs2_2(1,:,:)),2); e1=flip(e1);
e2 = mean(squeeze(fs2_2(2,:,:)),2); e2=flip(e2);

f1 = mean(squeeze(fs2_3(1,:,:)),2); f1=flip(f1);
f2 = mean(squeeze(fs2_3(2,:,:)),2); f2=flip(f2);

xx = [1,3,10,30,100,300,1000];

plot(xx,a1,'-ro')
hold on
plot(xx,d1,'-r*')

plot(xx,b1,'-go')
plot(xx,e1,'-g*')

plot(xx,c1,'-bo')
plot(xx,f1,'-b*')
