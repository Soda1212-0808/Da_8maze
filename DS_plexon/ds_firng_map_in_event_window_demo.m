% 示例数据
spike_times = [0.1, 0.5, 0.8, 1.2, 1.7, 2.0];  % 神经元放电时间点
time_windows = [0.0, 1.0; 1.0, 2.0];  % 时间窗，每行代表一个时间窗

% 获取时间窗的数量
num_windows = size(time_windows, 1);

% 使用逻辑索引创建掩码矩阵
is_in_window = (spike_times >= time_windows(:, 1)) & (spike_times < time_windows(:, 2));

% 使用矩阵乘法提取每个时间窗中的放电时间点
% 首先，将 is_in_window 转置，然后使用乘法
% 将 spike_times 重复 num_windows 次，然后元素相乘

spikes_in_windows = cell2mat(arrayfun(@(i) spike_times(is_in_window(i, :)), 1:num_windows, 'UniformOutput', false));

spike_indices = arrayfun(@(i) find(is_in_window(i, :)), 1:num_windows, 'UniformOutput', false);



% 显示结果
for i = 1:num_windows
    fprintf('时间窗 [%f, %f) 中的放电时间点: ', time_windows(i, 1), time_windows(i, 2));
    disp(spikes_in_windows{i});
end
