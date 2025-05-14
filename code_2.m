clc;
clear;
close all;
glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');
% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.03, 100, 0.001, 5);
imu = imuadderr(trj.imu, imuerr);
davp0 = avperrset([0.5;-0.5;20], 0.1, [1;1;3]);
ins = insinit(avpadderr(trj.avp0,davp0), ts);

