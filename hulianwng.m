%% 水下电磁波 vs 声波衰减对比仿真
clear; close all; clc;

% 参数设置
z = 0:0.1:100;          % 深度范围 0~100 米
E0 = 1;                 % 初始场强 (V/m)
sigma = 4;              % 海水电导率 (S/m)

% 电磁波频率设置
f_em = [1e6, 10e6];     % 1MHz 和 10MHz（典型水下无线电频率）
alpha_em = 0.0173 * sqrt(f_em * sigma); % 电磁波衰减系数 (dB/m)
alpha_em_np = alpha_em / 8.686;          % 转换为 Neper/m

% 声波参数（假设频率 10kHz，衰减系数远低于电磁波）
f_audio = 10e3;         % 10kHz 声波
alpha_audio = 0.001;    % 声波衰减系数 (dB/m) 
alpha_audio_np = alpha_audio / 8.686;

%% 计算场强衰减曲线
% 电磁波衰减
E_em1 = E0 * exp(-alpha_em_np(1) * z);   % 1MHz
E_em2 = E0 * exp(-alpha_em_np(2) * z);   % 10MHz

% 声波衰减（归一化到相同初始场强）
E_audio = E0 * exp(-alpha_audio_np * z);

%% 绘图：衰减曲线对比
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;

% 电磁波曲线
p1 = plot(z, E_em1, 'r-', 'LineWidth', 2, 'DisplayName', '电磁波 1MHz');
p2 = plot(z, E_em2, 'r--', 'LineWidth', 2, 'DisplayName', '电磁波 10MHz');

% 声波曲线
p3 = plot(z, E_audio, 'b-', 'LineWidth', 2, 'DisplayName', '声波 10kHz');

% 标注关键点
text(20, 0.8, sprintf('电磁波 10MHz\n衰减至 1%% 仅需 ~5m'), 'Color', 'r');
text(60, 0.6, sprintf('声波 10kHz\n100m 衰减 <10%%'), 'Color', 'b');

% 图形美化
xlabel('水下深度 (m)', 'FontSize', 12);
ylabel('归一化场强 E(z)/E_0', 'FontSize', 12);
title('水下电磁波 vs 声波衰减对比', 'FontSize', 14);
legend([p1, p2, p3], 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 11);
ylim([0, 1.1]);

%% 附加图：衰减系数 vs 频率
f = logspace(3, 9, 100);    % 频率范围：1kHz~1GHz
alpha_em_all = 0.0173 * sqrt(f * sigma);

figure('Color', 'w', 'Position', [100, 100, 800, 400]);
semilogx(f, alpha_em_all, 'r-', 'LineWidth', 2); hold on;
semilogx([f_audio, f_audio], [0, 0.1], 'b-', 'LineWidth', 2); % 声波参考线

% 标注
text(1e5, 50, '电磁波衰减系数', 'Color', 'r');
text(1e4, 0.005, '声波衰减系数 (~0.001 dB/m)', 'Color', 'b');
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('衰减系数 \alpha (dB/m)', 'FontSize', 12);
title('电磁波 vs 声波衰减系数随频率变化', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);