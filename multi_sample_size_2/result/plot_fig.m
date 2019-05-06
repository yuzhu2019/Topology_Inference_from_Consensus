clear;clc;
load('err1_36_01.mat');
err1_1 = err1;
load('err1_36_03.mat');
err1_2 = err1;
load('err1_36_05.mat');
err1_3 = err1;
load('err2_36_01.mat');
err2_1 = err2;
load('err2_36_03.mat');
err2_2 = err2;
load('err2_36_05.mat');
err2_3 = err2;

a1 = mean(squeeze(err1_1(1,:,:)),2); a1=flip(a1);
a2 = mean(squeeze(err1_1(2,:,:)),2); a2=flip(a2);

b1 = mean(squeeze(err1_2(1,:,:)),2); b1=flip(b1);
b2 = mean(squeeze(err1_2(2,:,:)),2); b2=flip(b2);

c1 = mean(squeeze(err1_3(1,:,:)),2); c1=flip(c1);
c2 = mean(squeeze(err1_3(2,:,:)),2); c2=flip(c2);

d1 = mean(squeeze(err2_1(1,:,:)),2); d1=flip(d1);
d2 = mean(squeeze(err2_1(2,:,:)),2); d2=flip(d2);

e1 = mean(squeeze(err2_2(1,:,:)),2); e1=flip(e1);
e2 = mean(squeeze(err2_2(2,:,:)),2); e2=flip(e2);

f1 = mean(squeeze(err2_3(1,:,:)),2); f1=flip(f1);
f2 = mean(squeeze(err2_3(2,:,:)),2); f2=flip(f2);

xx = [1,3,10,30,100,300,1000];

plot(xx,a1,'-ro')
hold on
plot(xx,d1,'-r*')

plot(xx,b1,'-go')
plot(xx,e1,'-g*')

plot(xx,c1,'-bo')
plot(xx,f1,'-b*')


