% 直接使用已有的带噪声数据
X_noisy = att_series(1:200,:)';  % 3×870，转置为 3×870（行：通道，列：时间）
X_noisy2 = att_series(201:400,:)';  % 3×870，转置为 3×870（行：通道，列：时间）

% LSTM 网络定义
layers = [ ...
    sequenceInputLayer(3)         % 输入3通道数据（姿态角）
    lstmLayer(50, 'OutputMode', 'sequence')  % LSTM 隐藏层
    fullyConnectedLayer(3)         % 输出3通道
    regressionLayer];

% 训练选项
options = trainingOptions('adam', ...
    'MaxEpochs', 50, ...
    'MiniBatchSize', 16, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', false, ...
    'Plots', 'training-progress');

% 使用原始数据作为目标进行训练
net = trainNetwork(X_noisy, X_noisy2, layers, options);

% 去噪
X_denoised = predict(net, X_noisy);

% 绘制原始 vs 去噪信号
figure;
for i = 1:3
    subplot(3,1,i);
    plot(X_noisy(i,:), 'r'); hold on;
    plot(X_denoised(i,:), 'b');
    legend('原始数据', '去噪数据');
    title(['Attitude Series Column ', num2str(i)]);
end
