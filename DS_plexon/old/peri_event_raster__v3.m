% 导入数据
% spike_times: 包含所有神经元放电时间的向量（单位：秒）
% event_times: 包含所有事件时间的向量（单位：秒）

spike_times = spike_whole; % 替换为实际文件名或数据
event_times = data.events{8, 1}.timestamps  ; % 替换为实际文件名或数据

% 设置事件窗口
pre_event_time = 5;  % 事件前1秒
post_event_time = 5; % 事件后1秒

% 初始化变量以存储每个事件窗口内的放电时间
all_event_spike_times = [];
trial_indices = [];
% 计算每个事件窗口内的放电时间
for i = 1:length(event_times)
    % 定义当前事件的时间窗口
    current_event_time = event_times(i);
    window_start = current_event_time - pre_event_time;
    window_end = current_event_time + post_event_time;
    
    % 找到在当前时间窗口内的放电时间
    spike_indices = find(spike_times >= window_start & spike_times <= window_end);
    aligned_spike_times = spike_times(spike_indices) - current_event_time;
    % 存储对齐后的放电时间
    all_event_spike_times = [all_event_spike_times ;aligned_spike_times];
    trial_indices = [trial_indices; i * ones(length(aligned_spike_times ), 1)];
end

% 设置直方图参数
bin_size = 0.05; % 直方图的bin大小（单位：秒）
edges = -pre_event_time:bin_size:post_event_time; % 直方图的边界

% 生成并绘制直方图
figure;
histogram(all_event_spike_times, edges);
xlabel('Time (s)');
ylabel('Spike Count');
title('Peri-Event Time Histogram');
grid on;


% 绘制 raster 图
figure;
scatter(all_event_spike_times, trial_indices, 10, 'k', 'filled');
xlabel('Time (s)');
ylabel('Trial');
title('Peri-Event Raster Plot');
xlim([-pre_event_time, post_event_time]);
ylim([0, length(event_times) + 1]);
grid on;

% 计算放电率
bin_counts = histcounts(all_event_spike_times, edges);
num_events = length(event_times);
firing_rate = bin_counts / (bin_size * num_events); % 归一化为放电率（单位：Hz）

% 绘制放电率直方图
figure;
bar(edges(1:end-1) + bin_size/2, firing_rate, 'histc');
xlabel('Time (s)');
ylabel('Firing Rate (Hz)');
title('Peri-Event Time Histogram with Firing Rate');
xlim([-pre_event_time, post_event_time]);
grid on;

