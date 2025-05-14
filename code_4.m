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

% KF filter
rk = poserrset([1;1;3]);% rk为卫导位置测量误差
kf = kfinit(ins, davp0, imuerr, rk);
% 设置最小协方差矩阵
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  kf.pconstrain=1;
len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);
timebar(nn, len, '15-state SINS/GPS Simulation.'); 
ki = 1;
% 进入滤波
for k=1:nn:len-nn+1
    k1 = k+nn-1;  
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);%定义状态转移矩阵
    kf = kfupdate(kf);
    if mod(t,1)==0 %ins间隔0.2s，gps间隔1s
        posGPS = trj.avp(k1,7:9)' + davp0(7:9).*randn(3,1);  % GPS pos simulation with some white noise
        kf = kfupdate(kf, ins.pos-posGPS, 'M'); 
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');%反馈补偿
        avp(ki,:) = [ins.avp', t];
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';  ki = ki+1;
    end
    timebar;
end
avp(ki:end,:) = [];  xkpk(ki:end,:) = []; 
% show results
insplot(avp);
avperr = avpcmpplot(trj.avp, avp);
kfplot(xkpk, avperr, imuerr);
