clear all
% 设置事件窗口

psth_window = [-5,5]
pre_event_time = 5;  % 事件前1秒
post_event_time = 5; % 事件后1秒
% 设置直方图参数
bin_size = 0.05; % 直方图的bin大小（单位：秒）


Path = 'H:\CA3_reprocess\each mice\';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};
for curr_animal=1:length(animals)
    % 获取文件夹中的所有内容
    contents = dir([Path animals{curr_animal}]);
    % 获取所有子文件夹的名称
    recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

    for curr_file=1:length(recording_files)
        %处理event数据
        event_file=dir(fullfile(Path, animals{curr_animal} ,recording_files{curr_file} ,'*data_m1.csv'));
        data_event=csvread(fullfile(event_file.folder , event_file.name));
        %         [event_idx,event_type]=findgroups(data_event(:,1));
        % 指定排序顺序
        desiredOrder = [26, 62, 17, 71, 35, 53];
        % 创建排序索引
        [~, sortIndex] = ismember(data_event(:, 1), desiredOrder);
        % 根据desiredOrder排序
        [~, sortOrder] = sort(sortIndex);
        % 对矩阵进行排序
        data_event_sorted = data_event(sortOrder, :);
        data_event_sorted_correct= data_event_sorted(data_event_sorted(:,2)==1,:);
        event_times=data_event_sorted_correct(:,4);


        neuron_files=dir(fullfile(Path ,animals{curr_animal} , recording_files{curr_file} ,'*.t64'));
        spike_whole = arrayfun(@(f) readmclusttfile(fullfile(f.folder, f.name))'/10000, ...
            neuron_files, 'UniformOutput', false);
        spike_times=spike_whole{1};

        %%整个spike的平均值和标准差
        all_counts = histcounts(spike_times, min(spike_times):bin_size:max(spike_times));
        mean_counts=mean(all_counts);
        sem_counts=std(all_counts);


        psth_bins = psth_window(1):bin_size:psth_window(2); % 直方图的边界


        % 使用 bsxfun 快速计算对齐的放电时间
        aligned_spike_times = bsxfun(@minus, spike_times, event_times');
        % 将对齐的放电时间限制在事件窗口内
        valid_spikes = aligned_spike_times(aligned_spike_times >= psth_window(1) & aligned_spike_times <= psth_window(2));
        % 为 raster 图准备数据
        [~,trial_indices ] = find(aligned_spike_times >= psth_window(1) & aligned_spike_times <= psth_window(2));
        %         % 计算每个 bin 内的放电次数
        %         bin_counts = histcounts(valid_spikes, psth_bins);
        %         bin_zscore=(bin_counts-mean_counts)/sem_counts;
        %         % 计算放电率
        %         firing_rate = bin_counts / (bin_size * length(event_times)); % 归一化为放电率（单位：Hz）
        stim_bins= event_times+ psth_bins;
        bin_counts = cell2mat(arrayfun(@(x) histcounts(spike_times, stim_bins(x,:)), (1:length(event_times))', 'UniformOutput', false));
        firing_rate=mean(bin_counts,1)/bin_size;
        firing_zscore=(mean(bin_counts,1)-mean_counts)/sem_counts;



    end
end


        %% 创建一个新图，放电率直方图
        figure;
        % 绘制 raster 图

        nexttile;
        colors_list = [0, 0, 0;    % 红色
            1, 0, 0;    % 蓝色
            0, 0, 1];   % 绿色

        % 初始化颜色数组
        colors = zeros(size(trial_indices, 1), 3);

        % 为每种颜色创建逻辑索引
        idx_{2} = ismember(data_event_sorted_correct(:,1), [26, 62]);
        idx_{1} = ismember(data_event_sorted_correct(:,1), [17, 71]);
        idx_{3} = ismember(data_event_sorted_correct(:,1), [35, 53]);

        % 找到匹配的trial_indices
        idx_trial_2 = ismember(trial_indices, find( idx_{2}));
        idx_trial_1 = ismember(trial_indices, find( idx_{1}));
        idx_trial_3 = ismember(trial_indices, find( idx_{3}));

        % 分配颜色
        colors(idx_trial_2, :) = repmat(colors_list(1, :), sum(idx_trial_2), 1);
        colors(idx_trial_1, :) = repmat(colors_list(2, :), sum(idx_trial_1), 1);
        colors(idx_trial_3, :) = repmat(colors_list(3, :), sum(idx_trial_3), 1);
        % 绘制散点图
        scatter(valid_spikes, trial_indices, 5, colors, 'filled');
        %         scatter(valid_spikes, trial_indices, 10, 'k', 'filled');
        xlabel('Time (s)');
        ylabel('Trial');
        title('Peri-Event Raster Plot');
        xlim(psth_window);
        ylim([0, length(event_times) + 1]);
     

        nexttile;
         hold on
        plot(psth_bins(1:end-1),smoothdata((mean(bin_counts(idx_{2},: ),1)-mean_counts)/sem_counts,'gaussian',20),'color', colors_list(1,:))
        plot(psth_bins(1:end-1),smoothdata((mean(bin_counts(idx_{1} ,:),1)-mean_counts)/sem_counts,'gaussian',20),'color', colors_list(2,:))
        plot(psth_bins(1:end-1),smoothdata((mean(bin_counts(idx_{3} ,:),1)-mean_counts)/sem_counts,'gaussian',20),'color', colors_list(3,:))

        % bar(edges(1:end-1) + bin_size/2, bin_zscore, 'histc');
        xlabel('Time (s)');
        %         ylabel('Firing Rate (Hz)');
        ylabel('zscore');
        title('Peri-Event Time Histogram with Firing Rate');
       xlim(psth_window);
 
