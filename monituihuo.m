T_init = 100; % 初始温度
T_min = 10; % 最小温度
alpha = 0.95; % 温度衰减系数
max_iter = 1000; % 最大迭代次数
current_n1 = 50000; % 初始 n_1
current_n2 = 20000; % 初始 n_2
current_angle = func_angle(current_n1, current_n2); % 初始目标函数值
adaptive_target = current_angle; % 初始时将目标设为当前的角度值
T = T_init;
best_n1 = current_n1;
best_n2 = current_n2;
best_angle = current_angle;

% 初始化适应性目标

iterationStore = zeros(1, max_iter); % 存储迭代次数
angleStore = zeros(1, max_iter); % 存储角度（误差）

% 迭代过程
for iter = 1:max_iter

    % 生成新的解
    new_n1 = current_n1 + round(100 * randn()); % 随机扰动
    new_n2 = current_n2 + round(100 * randn());
    % 检查 n_1 是否大于 n_2，如果不是则重新选择
    while new_n1 <= new_n2
        new_n1 = current_n1 + round(100 * randn());
        new_n2 = current_n2 + round(100 * randn());
    end
    % 计算新的目标函数值
    new_angle = func_angle(new_n1, new_n2);

    % 计算当前解与适应性目标之间的误差
    error_diff = (new_angle - adaptive_target)^2; % 均方误差

    % 计算目标函数的变化
    delta_angle = new_angle - current_angle;

    % 接受新解的条件：如果误差更小或根据退火准则接受
    if error_diff < (adaptive_target - current_angle)^2 || rand() < exp(-delta_angle / T)
        current_n1 = new_n1;
        current_n2 = new_n2;
        current_angle = new_angle;
    end

    % 更新最优解
    if current_angle < best_angle
        best_n1 = current_n1;
        best_n2 = current_n2;
        best_angle = current_angle;
    end

    % 更新适应性目标：通过误差平方和调整目标
    adaptive_target = 0.5 * adaptive_target + 0.5 * current_angle; % 基于历史优化更新目标

    % 降低温度
    T = T * alpha;

    % 存储每轮最优解
    iterationStore(iter) = iter;
    angleStore(iter) = best_angle;

    % 打印进度
    fprintf('Iteration %d, Temperature %.3f, Best ANGLE: %.5f, Adaptive Target: %.5f,angle each time:%5f\n', ...
        iter, T, best_angle, adaptive_target,new_angle);

    % 早停条件
    if T < T_min
        break;
    end
end

% 输出最优解
fprintf('Best n_1: %d, Best n_2: %d\n', best_n1, best_n2);






%% 原始算法
angleOrigin = zeros(20,1);
for i = 1:20
angleOrigin(i) = func_angle(60471,19614);
end
