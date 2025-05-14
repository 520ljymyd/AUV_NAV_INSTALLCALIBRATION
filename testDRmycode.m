%% Dead Recoding simulation with heading bias analysis
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS_DR.

clc; clear;
glvs;
trj = trjfile('trj10ms.mat'); % 获得 avp, imu, avp0, wat, ts, repeats
[nn, ts, nts] = nnts(2, trj.ts);

% 里程计安装参数
inst = [0;0;0.5]*glv.deg;  kod = 1;  qe = 0; dT = 0;
trjod = odsimu(trj, inst, kod, qe, dT, 0);
imuerr = imuerrset(0.01, 50, 0.001, 5);
imu = imuadderr(trjod.imu, imuerr);

davp = avperrset([0;0;0], 0, 0);  % 误差初始化

% 偏差参数
%dinst = [0;0;0]*glv.min; dkod = 0.25;
dinst = [0;0;0.5]*glv.deg; dkod = 0.05;

dr = drinit(avpadderr(trjod.avp0,davp), inst+dinst, kod*(1+dkod), ts);
len = length(imu); avp = prealloc(fix(len/nn), 10);
ki = timebar(nn, len, 'DR simulation.');

% 存储航向安装偏差角随时间变化
bias_angles = zeros(fix(len/nn), 1);
time_series = zeros(fix(len/nn), 1);

for k=1:nn:len-nn+1
    k1 = k+nn-1;
    wm = imu(k:k1,1:3);  dS= sum(trjod.od(k:k1,1)); t = imu(k1,end);
    dr = drupdate(dr, wm, dS);
    avp(ki,:) = [dr.avp; t]';
    
    % 计算当前时刻的航向安装偏差角
    position_currentDR = avp(ki, end-3:end-1);
    position_currentDR(1:2) = position_currentDR(1:2) * (180 / pi);
    
    position_currentTrue = trjod.avp(ki, end-3:end-1);
    position_currentTrue(1:2) = position_currentTrue(1:2) * (180 / pi);
    
    positionO = trjod.avp(1, end-3:end-1);
    positionO(1:2) = positionO(1:2) * (180 / pi);
    
    [~, angle_deg] = calculate_heading_bias(positionO(1:2), position_currentDR(1:2), position_currentTrue(1:2));
    
    % 记录数据
    bias_angles(ki) = angle_deg;
    time_series(ki) = t;
    
    ki = timebar;
end

dr.distance;

% 绘制安装偏差角随时间变化
% figure;
% plot(time_series, bias_angles, 'b-', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Heading Installation Bias Angle (deg)');
% title('Heading Bias Angle Convergence Analysis');
% grid on;
