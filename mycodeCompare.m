% SINS/GPS intergrated navigation simulation using kalman filter.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% aa2phi, aa2mu, aa2phimu
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011
glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');
% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.01, 1, 0.001, 1);
%imuerr = imuerrset(0.00, 0, 0.00, 0);
%imuerr = imuerrset(0.001,1);
imu = imuadderr(trj.imu, imuerr);
%davp0 = avperrset([0.5;-0.5;20], 0.1, [1;1;3]);
davp0 = avperrset([0;0;0], 0, [0;0;0]);
ins = insinit(avpadderr(trj.avp0,davp0), ts);
% KF filter
rk = poserrset([1;1;3]);
kf = kfinit(ins, davp0, imuerr, rk);
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  kf.pconstrain=1;
len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);
timebar(nn, len, '15-state SINS/GPS Simulation.');
ki = 1;
for k=1:nn:len-nn+1
    k1 = k+nn-1;
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);
    if mod(t,1)==0
        posGPS = trj.avp(k1,7:9)' + davp0(7:9).*randn(3,1);  % GPS pos simulation with some white noise
        %posGPS = trj.avp(k1,7:9)' ;  % GPS pos simulation with some white noise
        kf = kfupdate(kf, ins.pos-posGPS, 'M');
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        avp(ki,:) = [ins.avp', t];
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';  ki = ki+1;
    end
    timebar;
end
avp(ki:end,:) = [];  xkpk(ki:end,:) = [];
% show results
% insplot(avp);
% avperr = avpcmpplot(trj.avp, avp);
% kfplot(xkpk, avperr, imuerr);
%% DVL 获得d系速度，得到d系位移
inst=[0,0,0.5]*glv.deg;
kod = 1.005;
% trjod = odsimu(trj, inst,kod);
trjod = odsimu(trj, inst,kod);
v_od = trjod.avp(:, 4:6);
C_DVLb = a2mat(inst);
v_d = (C_DVLb * v_od')';  % 转置后矩阵乘法，获得DVL测量速度
% 噪声添加
sigma = 0.001;
gaussian_noise = sigma*randn(size(v_d));
v_d = v_d +gaussian_noise;
% 生成系统性偏置噪声
% bias_noise = 0.05 * randn(size(v_d));  % 偏置噪声
% v_d= v_d + bias_noise;
% shape = 2;   % 形状参数，控制噪声的“陡峭”程度
% scale = 0.1;  % 规模参数，控制噪声的大小
% gamma_noise = gamrnd(shape, scale, size(v_d));
% v_d = v_d + gamma_noise;

% v_d = v_d(1:100:end,:);
t_d = trjod.avp(:, end);  % 时间序列
% 初始化 d 系位移
s_d = zeros(size(v_d));
% 计算时间增量 dt_d
dt_d = diff(t_d);
dt_d = [dt_d; dt_d(end)];  % 确保时间序列长度一致
% 逐步累积计算 d 系位移
s_d = cumsum(v_d .* dt_d, 1);  % 逐行累积计算位移



%% 获取v_b以及s_b
v_n = avp(:, 4:6);% 提取导航系速度 (v_n)
% 获取姿态角，并计算方向余弦矩阵 (Cnb)
att1 = avp(:, 1:3);
Cnb = a2mat(att1);
% 计算载体坐标系速度 (v_b)
v_b = (Cnb'* v_n')';  % 转置矩阵后进行乘法，得到载体系速度
% 初始化b系位移
s_b = zeros(size(v_b));
% 时间间隔增量dt
dt_b= diff(avp(:, end));  % 计算时间增量
dt_b = [dt_b; dt_b(end)];  % 确保时间序列长度一致
% 遍历每一时刻，计算 b 系位移

s_b = cumsum(v_b .* dt_b, 1);  % 逐行累积计算位移

% % 创建新的时间序列  线性插值
t_b = avp(:,end);
t_b_new = linspace(t_b(1), t_b(end), length(t_d));

% 90000个数据点的时间序列

% 使用插值函数将s_b
s_b = interp1(t_b, s_b, t_b_new, 'linear');

%
%     figure;
%     hold on;
%     plot(t_d,s_d(:,1),'r',LineWidth=0.5);
%     plot(t_d,s_d(:,2),'b',LineWidth=0.5);
%     plot(t_d,s_d(:,3),'g',LineWidth=0.5);
%     xlabel('Time (s)');
%     ylabel('Displacement (m)');
%     title('Displacement in x,y,z direction   (D)');
%
% figure;
% hold on;
% plot(t_b_new,s_b(:,1),'r',LineWidth=1.5);
% plot(t_b_new,s_b(:,2),'b',LineWidth=1.5);
% plot(t_b_new,s_b(:,3),'g',LineWidth=1.5);
% xlabel('Time (s)');
% ylabel('Displacement (m)');
% title('Displacement in x,y,z direction    (B)');
% 计算 b 系 和 d 系的位移矩阵
Sb = s_b';  % 3×N 矩阵 (b 系位移)
% Vb = v_b';
Sd = s_d'; % 3×N 矩阵 (d 系位移)
% Vd = v_d';


% 初始化角度记录数组
anglesOrigin = zeros(length(s_b), 1);  % 创建一个与时间序列长度相同的数组用于存储角度

% 遍历每一时刻，计算安装偏差角度并记录
for k = 2:length(s_b)
    % 计算安装偏差矩阵 C_b^d
    Sb2 = Sb(:, k);  % Sb 的最后一列
    %   Sb1 = sum(Sb(:,2:k), 2);
    %   Sb2 = Sb(:,end-5000);  % Sb 的倒数第二列
    Sb1 = sum(Sb(:,1:k-1), 2);

    Sd2 = Sd(:, k);  % Sb 的最后一列
    %Sd2 = Sd(:, end-5000);  % Sb 的倒数第二列
    Sd1 = sum(Sd(:,1:k-1), 2);

    % 计算 Sb1 和 Sb2 的叉积
    cross_result = cross(Sb1, Sb2);
    normalized_cross = cross_result / norm(cross_result);
    cross_result2 = cross(cross_result, Sb1);
    normalized_cross2 = cross_result2 / norm(cross_result2);

    % 创建 v1 向量
    v1 = [Sb1/norm(Sb1), normalized_cross, normalized_cross2];

    % 计算 Sd1 和 Sd2 的叉积
    cross_result_sd = cross(Sd1, Sd2);
    normalized_cross_sd = cross_result_sd / norm(cross_result_sd);
    cross_result_sd2 = cross(cross_result_sd, Sd1);
    normalized_cross_sd2 = cross_result_sd2 / norm(cross_result_sd2);

    % 创建 v2 向量
    v2 = [Sd1/norm(Sd1), normalized_cross_sd, normalized_cross_sd2];

    % 计算安装偏差矩阵 C_b^d
    Cbd = v1 * v2';  % 计算 v1 乘以 v2 的转置

    % 计算偏差角度
    angle = Cbd(2,1) * 180 / pi;

    % 记录角度
    anglesOrigin(k) = angle;
end

figure;
t = 1:length(anglesOrigin);
plot(t, anglesOrigin, 'b', 'LineWidth', 1.5);  % Plot the first set of angles (blue line)
xlabel('Time (s)');
ylabel('Installation Bias Angle (degrees)');
title('Installation Bias Angle vs Time');
grid on;
yline(0.5, '--r', 'LineWidth', 2);
ylim([-0.5,0.7]);
legend('Algorithm 1', 'Algorithm 2');






