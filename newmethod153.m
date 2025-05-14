%这个代码是新标定算法，新标定算法
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
len = length(imu); [avp, xkpk] = prealloc(fix(len/1), 10, 2*kf.n+1);
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
inst=[0,0,0.2]*glv.deg;
% kod = 1.005;
trjod = odsimu(trj, inst);
% trjod = odsimu(trj, inst);
v_od = trjod.avp(:, 4:6);
C_DVLb = a2mat(inst);
v_d = (C_DVLb * v_od')';  % 转置后矩阵乘法，获得DVL测量速度
% 噪声添加
% sigma = 0.001;
% gaussian_noise = sigma*randn(size(v_d));
% v_d = v_d +gaussian_noise;
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


% figure;
% hold on;
% plot(t_d,s_d(:,1),'r',LineWidth=0.5);
% plot(t_d,s_d(:,2),'b',LineWidth=0.5);
% plot(t_d,s_d(:,3),'g',LineWidth=0.5);
% xlabel('Time (s)');
% ylabel('Displacement (m)');
% title('Displacement in x,y,z direction   (D)');

% figure;
% hold on;
% plot(t_b_new,s_b(:,1),'r',LineWidth=1.5);
% plot(t_b_new,s_b(:,2),'b',LineWidth=1.5);
% plot(t_b_new,s_b(:,3),'g',LineWidth=1.5);
% xlabel('Time (s)');
% ylabel('Displacement (m)');
% title('Displacement in x,y,z direction    (B)');
%% 计算 b 系 和 d 系的位移矩阵
Sb = s_b';  % 3×N 矩阵 (b 系位移)
% Vb = v_b';
Sd = s_d'; % 3×N 矩阵 (d 系位移)
% Vd = v_d';

%%
% 初始化角度记录数组
anglesOrigin = zeros(length(s_b), 1);  % 创建一个与时间序列长度相同的数组用于存储角度

% 遍历每一时刻，计算安装偏差角度并记录
for k = 10000:length(s_b)
    % 计算安装偏差矩阵 C_b^d

    Sb1 = Sb(:, k);  % Sb 的最后一列
    Sb2 = sum(Sb(:,k-n_1:k-n_2), 2);
    Sd1 = Sd(:, k);  % Sb 的最后一列
    Sd2 = sum(Sd(:,k-n_1:k-n_2), 2);

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
plot(t, anglesOrigin, 'b', 'LineWidth', 0.7);  % Plot the first set of angles (blue line)
xlabel('Time (s)');
ylabel('Installation Bias Angle (degrees)');
title('Installation Bias Angle vs Time');
grid on;
yline(inst(3)*180/pi, '--r', 'LineWidth', 2);
ylim([-0.5,0.7]);
legend('Algorithm 1', 'Algorithm 2');

%% 智能算法   模拟退火算法
%angle = func_angle(n_1,n_2);
% 定义粒子群的参数
% nvars = 2;  % 我们有两个优化参数 n_1 和 n_2
% lb = [500,100];  % n_1 和 n_2 的下界（可以根据实际情况调整）
% ub = [1000,500];  % n_1 和 n_2 的上界（可以根据实际情况调整）
%
% % 设置优化选项
% options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100, 'Display', 'iter');
%
% % 运行粒子群优化
% [best_params, best_error] = particleswarm(@optimization_function, nvars, lb, ub, options);
%
% % 输出最优解
% n_1_opt = round(best_params(1));  % 最优的 n_1
% n_2_opt = round(best_params(2));  % 最优的 n_2
% fprintf('最优参数: n_1 = %d, n_2 = %d\n', n_1_opt, n_2_opt);
% fprintf('最小误差: %f\n', best_error);
% 初始化图形窗口

% T_init = 100;
% T_min = 10;
% alpha = 0.95;
% max_iter = 1000; % 最大迭代次数
% current_n1 = 50000; % 初始 n_1
% current_n2 = 1000; % 初始 n_2
% current_angle = func_angle(current_n1,current_n2);
% T = T_init;
% best_n1 = current_n1;
% best_n2 = current_n2;
% best_angle = current_angle;
%
% iterationStore = zeros(1, max_iter); % Store iterations
% angleStore = zeros(1, max_iter);
% for iter = 1:max_iter
%     % 生成新的解
%     new_n1 = current_n1 + round(100 * randn()); % 随机扰动
%     new_n2 = current_n2 + round(100 * randn());
%
%     % 目标函数值
%     new_angle = func_angle(new_n1, new_n2);
%
%     % 计算目标函数的变化
%     delta_angle = new_angle - current_angle;
%
%     % 接受新解的条件
%     if delta_angle < 0 || rand() < exp(-delta_angle / T)
%         current_n1 = new_n1;
%         current_n2 = new_n2;
%         current_angle = new_angle;
%     end
%
%     % 更新最优解
%     if abs(current_angle-optnum) < abs(best_angle-optnum)
%         best_n1 = current_n1;
%         best_n2 = current_n2;
%         best_angle = current_angle;
%     end
%
%     % 降低温度
%     T = T * alpha;
%     iterationStore(iter) = iter;
%     angleStore(iter) = best_angle;
%
%     % 打印进度
%     fprintf('Iteration %d, Temperature %.3f, Best ANGLE: %.5f\n', iter, T, best_angle);
%
%     % 早停条件
%     if T < T_min
%         break;
%     end
% end
%
% % 输出最优解
% fprintf('Best n_1: %d, Best n_2: %d\n', best_n1, best_n2);



