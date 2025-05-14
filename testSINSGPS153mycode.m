% 初始化
glvs;
psinstypedef(153);
trj = trjfile('trj10ms.mat');
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.01, 1, 0.001, 1);
imu = imuadderr(trj.imu, imuerr);
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
        kf = kfupdate(kf, ins.pos-posGPS, 'M');
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        avp(ki,:) = [ins.avp', t];
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';  ki = ki+1;
    end
    timebar;
end
avp(ki:end,:) = [];  xkpk(ki:end,:) = [];

% 获取v_b以及s_b
v_n = avp(:, 4:6); % 提取导航系速度 (v_n)
att1 = avp(:, 1:3);
Cnb = a2mat(att1);
v_b = (Cnb'* v_n')';  % 转置矩阵后进行乘法，得到载体系速度
sigma = 0.1;
gaussian_noise = sigma*randn(size(v_b));
v_b = v_b + gaussian_noise;

% 初始化b系位移
s_b = zeros(size(v_b));
dt_b = diff(avp(:, end));  % 计算时间增量
dt_b = [dt_b; dt_b(end)];  % 确保时间序列长度一致
for k = 2:length(v_b)
    s_b(k, :) = s_b(k-1, :) + v_b(k, :) * dt_b(k);  % 累积位移
end

% 创建新的时间序列，90000个数据点，线性插值
t_b = avp(:,end);
t_b_new = linspace(t_b(1), t_b(end), length(imu));  % 90000个数据点的时间序列
s_b = interp1(t_b, s_b, t_b_new, 'linear');
v_b = interp1(t_b, v_b, t_b_new, 'linear');

%% DVL 获得d系速度，得到d系位移
inst = [0.0, -3.3, 2.2] * glv.deg;
kod = 1.5;
trjod = odsimu(trj, inst, kod);
v_od = trjod.avp(:, 4:6);
C_DVLb = a2mat(inst);
v_d = (C_DVLb * v_od')';  % 转置后矩阵乘法，获得DVL测量速度
sigma = 0.1;
gaussian_noise = sigma*randn(size(v_d));
v_d = v_d + gaussian_noise;
t_d = trjod.avp(:, end);  % 时间序列
s_d = zeros(size(v_d));
dt_d = diff(t_d);
dt_d = [dt_d; dt_d(end)];  % 确保时间序列长度一致
for k = 2:length(v_d)
    s_d(k, :) = s_d(k-1, :) + v_d(k, :) * dt_d(k);
end

%% 卡尔曼滤波迭代估计DVL安装角
% 初始化卡尔曼滤波参数
X = zeros(3, 1);  % 状态变量 [α, β, γ]
P = eye(3) * 1000;  % 初始协方差矩阵
Q = eye(3) * 0.001;  % 过程噪声协方差
R = eye(3) * 0.01;  % 观测噪声协方差
C = eye(3);  % 初始转换矩阵

% 迭代估计
max_iter = 10;
tol = 0.01;
for iter = 1:max_iter
    for k = 1:length(v_b)
        % 计算观测残差 Z = C * v_d - v_b
        Z = C * v_d(k, :)' - v_b(k, :)';
        
        % 观测矩阵 H
        H = [
            0,      v_d(k, 3), -v_d(k, 2);
            -v_d(k, 3), 0,      v_d(k, 1);
            v_d(k, 2), -v_d(k, 1), 0
        ];
        
        % 卡尔曼增益计算
        S = H * P * H' + R;
        K = P * H' / S;
        
        % 状态更新
        X = X + K * Z;
        P = (eye(3) - K * H) * P;
        
        % 过程噪声更新
        P = P + Q;
        
        % 更新转换矩阵
        alpha = X(1); beta = X(2); gama = X(3);
        C = [
            1,   gama,  -beta;
            -gama,  1,   alpha;
            beta,  -alpha,   1
        ];
    end
    
    % 检查收敛条件
    if norm(X) < tol
        break;
    end
end

% 输出估计结果
fprintf('Estimated installation angles (degrees):\n');
fprintf('α: %.4f, β: %.4f, γ: %.4f\n', rad2deg(X(1)), rad2deg(X(2)), rad2deg(X(3)));

% 验证结果
v_b_estimated = (C * v_d')';
errors = v_b_estimated - v_b;
error_norms = sqrt(sum(errors.^2, 2));

% 绘制误差曲线
figure;
plot(error_norms);
xlabel('Sample Index');
ylabel('Error Norm (m/s)');
title('DVL Calibration Validation');
grid on;

% 打印统计结果
fprintf('Max Error: %.4f m/s\n', max(error_norms));
fprintf('Mean Error: %.4f m/s\n', mean(error_norms));
fprintf('Std Error: %.4f m/s\n', std(error_norms));