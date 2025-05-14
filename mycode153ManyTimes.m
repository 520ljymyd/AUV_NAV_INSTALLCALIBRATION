glvs;
angles = [];  % 用于存储每次计算的角度
% 使用0.5秒的时间步长循环
for Time = 0:50:900
    ts = 0.01;  % 采样间隔为0.01（100Hz）

    % 初始avp
    avp0 = [[0;0;-30]; [0;0;0]; [45.7796*glv.deg; 126.6705*glv.deg; -50]];  % 初始化avp
    % 设置轨迹段
    xxx = [];
    seg = trjsegment(xxx, 'init', 0);  % 初始化
    seg = trjsegment(seg, 'accelerate', 3, xxx, 1);  % 加速
    seg = trjsegment(seg, 'uniform', Time - 1);  % 均匀轨迹
    seg = trjsegment(seg, 'coturnleft', 45, 2, xxx, 4);  % 左转
    seg = trjsegment(seg, 'uniform', 100);  % 继续均匀轨迹
    trj = trjsimu(avp0, seg.wat, ts, 1);
   % trjfile('trj10ms.mat', trj);
    close all;
    % 使用kalman滤波更新轨迹
    psinstypedef(153);
    %trj = trjfile('trj10ms.mat');
    % 初始设置
    [nn, ts, nts] = nnts(2, trj.ts);
    imuerr = imuerrset(0.01, 1, 0.001, 1);
    imu = imuadderr(trj.imu, imuerr);
    davp0 = avperrset([0; 0; 0], 0, [0; 0; 0]);
    ins = insinit(avpadderr(trj.avp0, davp0), ts);

    % KF滤波器设置
    rk = poserrset([1; 1; 3]);
    kf = kfinit(ins, davp0, imuerr, rk);
    kf.Pmin = [avperrset(0.01, 1e-4, 0.1); gabias(1e-3, [1, 10])].^2;
    kf.pconstrain = 1;
    len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);
%     timebar(nn, len, '15-state SINS/GPS Simulation.');
    ki = 1;

    % 循环处理每个时间点的IMU数据
    for k = 1:nn:len - nn + 1
        k1 = k + nn - 1;
        wvm = imu(k:k1, 1:6);  t = imu(k1, end);
        ins = insupdate(ins, wvm);
        kf.Phikk_1 = kffk(ins);
        kf = kfupdate(kf);

        if mod(t, 1) == 0
            posGPS = trj.avp(k1, 7:9)' + davp0(7:9) .* randn(3, 1);  % GPS位置模拟带噪声
            kf = kfupdate(kf, ins.pos - posGPS, 'M');
            [kf, ins] = kffeedback(kf, ins, 1, 'avp');
            avp(ki, :) = [ins.avp', t];
            xkpk(ki, :) = [kf.xk; diag(kf.Pxk); t]';
            ki = ki + 1;
        end
%         timebar;
    end
    avp(ki:end, :) = [];  xkpk(ki:end, :) = [];
    close all;
    % 获取b系速度和位移
    v_n = avp(:, 4:6);  % 提取导航系速度(v_n)
    att1 = avp(:, 1:3);
    Cnb = a2mat(att1);
    v_b = (Cnb' * v_n')';  % 转置矩阵后进行乘法得到载体系速度
    s_b = zeros(size(v_b));  % 初始化b系位移
    dt_b = diff(avp(:, end));  % 计算时间增量
    dt_b = [dt_b; dt_b(end)];  % 确保时间序列一致
    for k = 2:length(v_b)
        s_b(k, :) = s_b(k-1, :) + v_b(k, :) * dt_b(k);  % 累积位移
    end
    t_b = avp(:, end);
    t_b_new = linspace(t_b(1), t_b(end), 90200);
    s_b = interp1(t_b, s_b, t_b_new, 'linear');

    % DVL 获得d系速度并得到d系位移
    inst = [0, 0, 0.5] * glv.deg;
    kod = 1.6;
    trjod = odsimu(trj, inst, kod);
    v_od = trjod.avp(:, 4:6);
    C_DVLb = a2mat(inst);
    v_d = (C_DVLb * v_od')';  % 转置后矩阵乘法获得DVL速度
    t_d = trjod.avp(:, end);
    s_d = zeros(size(v_d));  % 初始化d系位移
    dt_d = diff(t_d);
    dt_d = [dt_d; dt_d(end)];  % 确保时间序列一致
    for k = 2:length(v_d)
        s_d(k, :) = s_d(k-1, :) + v_d(k, :) * dt_d(k);  % 累积位移
    end

    % 计算安装偏差矩阵C_b^d
    Sb = s_b';
    Sd = s_d';

    Sb1 = Sb(:, end);  % Sb的最后一列
    Sb2 = sum(Sb(:, 1:end-1), 2);  % Sb的其他列的和
    Sd1 = Sd(:, end);  % Sd的最后一列
    Sd2 = sum(Sd(:, 1:end-1), 2);  % Sd的其他列的和

    % 计算Sb1和Sb2的叉积
    cross_result = cross(Sb1, Sb2);
    normalized_cross = cross_result / norm(cross_result);
    cross_result2 = cross(cross_result, Sb1);
    normalized_cross2 = cross_result2 / norm(cross_result2);
    v1 = [Sb1/norm(Sb1), normalized_cross, normalized_cross2];

    % 计算Sd1和Sd2的叉积
    cross_result_sd = cross(Sd1, Sd2);
    normalized_cross_sd = cross_result_sd / norm(cross_result_sd);
    cross_result_sd2 = cross(cross_result_sd, Sd1);
    normalized_cross_sd2 = cross_result_sd2 / norm(cross_result_sd2);
    v2 = [Sd1/norm(Sd1), normalized_cross_sd, normalized_cross_sd2];
    glvs;
    angles = [];  % 用于存储每次计算的角度
    % 使用0.5秒的时间步长循环
    for Time = 1:10:300
        ts = 0.01;  % 采样间隔为0.01（100Hz）

        % 初始avp
        avp0 = [[0;0;-30]; [0;0;0]; [45.7796*glv.deg; 126.6705*glv.deg; -50]];  % 初始化avp

        % 设置轨迹段
        xxx = [];
        seg = trjsegment(xxx, 'init', 0);  % 初始化
        seg = trjsegment(seg, 'accelerate', 3, xxx, 1);  % 加速
        seg = trjsegment(seg, 'uniform', Time - 1);  % 均匀轨迹
        seg = trjsegment(seg, 'coturnleft', 45, 2, xxx, 4);  % 左转
        seg = trjsegment(seg, 'uniform', 100);  % 继续均匀轨迹
        trj = trjsimu(avp0, seg.wat, ts, 1);
        trjfile('trj10ms.mat', trj);
        close all;
        % 使用kalman滤波更新轨迹
        psinstypedef(153);
        trj = trjfile('trj10ms.mat');

        % 初始设置
        [nn, ts, nts] = nnts(2, trj.ts);
        imuerr = imuerrset(0.01, 1, 0.001, 1);
        imu = imuadderr(trj.imu, imuerr);
        davp0 = avperrset([0; 0; 0], 0, [0; 0; 0]);
        ins = insinit(avpadderr(trj.avp0, davp0), ts);

        % KF滤波器设置
        rk = poserrset([1; 1; 3]);
        kf = kfinit(ins, davp0, imuerr, rk);
        kf.Pmin = [avperrset(0.01, 1e-4, 0.1); gabias(1e-3, [1, 10])].^2;
        kf.pconstrain = 1;
        len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);
        timebar(nn, len, '15-state SINS/GPS Simulation.');
        ki = 1;

        % 循环处理每个时间点的IMU数据
        for k = 1:nn:len - nn + 1
            k1 = k + nn - 1;
            wvm = imu(k:k1, 1:6);  t = imu(k1, end);
            ins = insupdate(ins, wvm);
            kf.Phikk_1 = kffk(ins);
            kf = kfupdate(kf);

            if mod(t, 1) == 0
                posGPS = trj.avp(k1, 7:9)' + davp0(7:9) .* randn(3, 1);  % GPS位置模拟带噪声
                kf = kfupdate(kf, ins.pos - posGPS, 'M');
                [kf, ins] = kffeedback(kf, ins, 1, 'avp');
                avp(ki, :) = [ins.avp', t];
                xkpk(ki, :) = [kf.xk; diag(kf.Pxk); t]';
                ki = ki + 1;
            end
            timebar;
        end
        avp(ki:end, :) = [];  xkpk(ki:end, :) = [];
        close all;
        % 获取b系速度和位移
        v_n = avp(:, 4:6);  % 提取导航系速度(v_n)
        att1 = avp(:, 1:3);
        Cnb = a2mat(att1);
        v_b = (Cnb' * v_n')';  % 转置矩阵后进行乘法得到载体系速度
        s_b = zeros(size(v_b));  % 初始化b系位移
        dt_b = diff(avp(:, end));  % 计算时间增量
        dt_b = [dt_b; dt_b(end)];  % 确保时间序列一致
        for k = 2:length(v_b)
            s_b(k, :) = s_b(k-1, :) + v_b(k, :) * dt_b(k);  % 累积位移
        end
        t_b = avp(:, end);
        t_b_new = linspace(t_b(1), t_b(end), 90200);
        s_b = interp1(t_b, s_b, t_b_new, 'linear');

        % DVL 获得d系速度并得到d系位移
        inst = [0, 0, 0.5] * glv.deg;
        kod = 1.6;
        trjod = odsimu(trj, inst, kod);
        v_od = trjod.avp(:, 4:6);
        C_DVLb = a2mat(inst);
        v_d = (C_DVLb * v_od')';  % 转置后矩阵乘法获得DVL速度
        t_d = trjod.avp(:, end);
        s_d = zeros(size(v_d));  % 初始化d系位移
        dt_d = diff(t_d);
        dt_d = [dt_d; dt_d(end)];  % 确保时间序列一致
        for k = 2:length(v_d)
            s_d(k, :) = s_d(k-1, :) + v_d(k, :) * dt_d(k);  % 累积位移
        end

        % 计算安装偏差矩阵C_b^d
        Sb = s_b';
        Sd = s_d';

        Sb1 = Sb(:, end);  % Sb的最后一列
        Sb2 = sum(Sb(:, 1:end-1), 2);  % Sb的其他列的和
        Sd1 = Sd(:, end);  % Sd的最后一列
        Sd2 = sum(Sd(:, 1:end-1), 2);  % Sd的其他列的和

        % 计算Sb1和Sb2的叉积
        cross_result = cross(Sb1, Sb2);
        normalized_cross = cross_result / norm(cross_result);
        cross_result2 = cross(cross_result, Sb1);
        normalized_cross2 = cross_result2 / norm(cross_result2);
        v1 = [Sb1/norm(Sb1), normalized_cross, normalized_cross2];

        % 计算Sd1和Sd2的叉积
        cross_result_sd = cross(Sd1, Sd2);
        normalized_cross_sd = cross_result_sd / norm(cross_result_sd);
        cross_result_sd2 = cross(cross_result_sd, Sd1);
        normalized_cross_sd2 = cross_result_sd2 / norm(cross_result_sd2);
        v2 = [Sd1/norm(Sd1), normalized_cross_sd, normalized_cross_sd2];

        % 计算C_b^d并得出偏差角
        Cbd = v1 * v2';  % 计算v1和v2的转置
        angle = Cbd(2, 1) * 180 / pi;  % 将结果转换为度

        % 存储计算结果
        angles = [angles; Time, angle];
        close all;
    end

    % 绘制偏差角随时间变化的曲线
    figure;
    plot(angles(:, 1), angles(:, 2));
    xlabel('时间 (秒)');
    ylabel('偏差角 (度)');
    title('安装偏差角随时间变化');
    grid on;

    % 计算C_b^d并得出偏差角
    Cbd = v1 * v2';  % 计算v1和v2的转置
    angle = Cbd(2, 1) * 180 / pi;  % 将结果转换为度

    % 存储计算结果
    angles = angle;
    close all;
end

% 绘制偏差角随时间变化的曲线
% figure;
% plot(angles(:, 1), angles(:, 2));
% xlabel('时间 (秒)');
% ylabel('偏差角 (度)');
% title('安装偏差角随时间变化');
% grid on;
