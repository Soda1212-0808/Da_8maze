% 创建样例数据
n = 1000; % 数据点的数量
sampling_frequency = 30; % 采样频率为10Hz
time_interval = 1 / sampling_frequency; % 每个数据点之间的时间间隔为0.1秒
x_data = linspace(0, 2*pi, n); % 在 0 到 2*pi 之间生成 n 个点
y_data = sin(x_data); % 计算正弦值
trajectory = [x_data; y_data]; % 生成 2*n 数组

% 提取 X 和 Y 坐标
x = trajectory(1, :);
y = trajectory(2, :);

% 创建一个新图窗
figure;
hold on; % 保持图窗中的图形
xlabel('X 坐标');
ylabel('Y 坐标');
title('点的轨迹逐渐绘制');
grid on; % 打开网格
axis equal; % 设置坐标轴比例相同
xlim([min(x)-1, max(x)+1]); % 设置 X 轴范围
ylim([min(y)-1, max(y)+1]); % 设置 Y 轴范围

% 动态绘制轨迹
for i = 1:length(x)
    if i == 1
        % 初始点
        h = plot(x(i), y(i), 'bo'); % 绘制第一个点
    else
        % 更新轨迹
        set(h, 'XData', x(1:i), 'YData', y(1:i)); % 更新已有的点
        plot(x(i), y(i), 'bo'); % 绘制新点
    end
    pause(time_interval / 100); % 按100倍速度暂停，以创建动画效果
end

hold off;
