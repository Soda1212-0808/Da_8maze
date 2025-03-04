clear all
% 设置事件窗口

psth_window = [-5,5];

% 设置直方图参数
bin_size = 0.05; % 直方图的bin大小（单位：秒）
smooth_window=10;
Path='E:\SDdata\WT\each mice\'
% Path = 'H:\CA3_reprocess\each mice\';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};
newfolderName = 'PSTH';
if exist(fullfile(Path,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end
for curr_animal=1:length(animals)
    % 获取文件夹中的所有内容
    animal=animals{curr_animal};
    contents = dir([Path animal]);
    % 获取所有子文件夹的名称
    recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

    for curr_file=1:2
        %处理event数据
        event_file=dir(fullfile(Path, animal ,recording_files{curr_file} ,'*data_m1.csv'));
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


        event_times{1}=data_event_sorted(:,3);
        event_times{2}=data_event_sorted(:,4);
        event_times{3}=data_event_sorted(:,6);
        event_times{4}=data_event_sorted(:,7);
        event_times{5}=data_event_sorted(:,8);
        event_times{6}=data_event_sorted(:,10);


        neuron_files=dir(fullfile(Path ,animal , recording_files{curr_file} ,'*.t64'));
        spike_whole = arrayfun(@(f) readmclusttfile(fullfile(f.folder, f.name))'/10000000000000, ...
            neuron_files, 'UniformOutput', false);

        for curr_cell=1:length(spike_whole)
            spike_times=spike_whole{curr_cell};
            spike_name=neuron_files(curr_cell).name(1:end-4);
            %%整个spike的平均值和标准差
            all_counts = histcounts(spike_times, min(spike_times):bin_size:max(spike_times));
            mean_counts=mean(all_counts);
            sem_counts=std(all_counts);
            psth_bins = psth_window(1):bin_size:psth_window(2); % 直方图的边界



            stim_bins=cell(1,length(event_times));
            trial_indices=cell(1,length(event_times));
            bin_counts=cell(1,length(event_times));
            valid_spikes=cell(1,length(event_times));
            for curr_event=1:length(event_times)
                % 使用 bsxfun 快速计算对齐的放电时间
                aligned_spike_times = bsxfun(@minus, spike_times, event_times{curr_event}');
                % 将对齐的放电时间限制在事件窗口内
                valid_spikes{curr_event} = aligned_spike_times(aligned_spike_times >= psth_window(1) & aligned_spike_times <= psth_window(2));
                % 为 raster 图准备数据
                [~,trial_indices{curr_event} ] = find(aligned_spike_times >= psth_window(1) & aligned_spike_times <= psth_window(2));
                %         % 计算每个 bin 内的放电次数
                %         bin_counts = histcounts(valid_spikes, psth_bins);
                %         bin_zscore=(bin_counts-mean_counts)/sem_counts;
                %         % 计算放电率
                %         firing_rate = bin_counts / (bin_size * length(event_times)); % 归一化为放电率（单位：Hz）
                stim_bins{curr_event}= event_times{curr_event}+ psth_bins;
                bin_counts{curr_event} = cell2mat(arrayfun(@(x) histcounts(spike_times, stim_bins{curr_event}(x,:)), (1:length(event_times{curr_event}))', 'UniformOutput', false));
                %             firing_rate=mean(bin_counts,1)/bin_size;
                %             firing_zscore=(mean(bin_counts,1)-mean_counts)/sem_counts;
            end

            %创建一个新图，放电率直方图
            figure('Position',[50 100 length(event_times)*300 800]);
            tt = tiledlayout(1,length(event_times),'TileSpacing','tight');
            title(tt,[ animal '-day-' num2str(curr_file) ' ' spike_name ])
            figure_name={'sample begin','sample reward','delay','choice begin','sample reward','end'}
            % 绘制 raster 图
            for curr_figure=1:length(event_times)
                t_figure = tiledlayout(tt,2,1);
                t_figure.Layout.Tile = curr_figure;
                title(t_figure,figure_name{curr_figure});

                nexttile(t_figure);
                colors_list = [0, 0, 0; 0.5, 0.5, 0.5;   % 红色
                    1, 0, 0;  1, 0.5, 0.5;    % 蓝色
                    0, 0, 1;0.5, 0.5, 1];   % 绿色

                % 初始化颜色数组
                colors = zeros(size(trial_indices{curr_figure}, 1), 3);
                event_name=[26 62 17 71 35 53];
                % 为每种event创建逻辑索引
                for curr_arm=1:6
                    idx_{curr_arm} = ismember(data_event_sorted(:,1), event_name(curr_arm));
                    % 找到匹配的trial_indices
                    idx_trial{curr_arm} = ismember(trial_indices{curr_figure}, find( idx_{curr_arm} ));
                    % 分配颜色
                    colors(idx_trial{curr_arm}, :) = repmat(colors_list(curr_arm, :), sum(idx_trial{curr_arm}), 1);
                end

                % 绘制散点图
                scatter(valid_spikes{curr_figure}, trial_indices{curr_figure}, 3, colors, 'filled');
                %         scatter(valid_spikes, trial_indices, 10, 'k', 'filled');
                xlabel('Time (s)');
                ylabel('Trial');
                title('Peri-Event Raster Plot');
                xlim(psth_window);
                ylim([0, length(event_times{curr_figure}) + 1]);


                nexttile(t_figure);

                hold on
                arrayfun(@(x) plot(psth_bins(1:end-1),smoothdata((mean(bin_counts{curr_figure}(idx_{x},: ),1)-mean_counts)/sem_counts,'gaussian',smooth_window),'LineWidth',1.5,'color', colors_list(x,:)),1:6)

                % bar(edges(1:end-1) + bin_size/2, bin_zscore, 'histc');
                xlabel('Time (s)');
                %         ylabel('Firing Rate (Hz)');
                ylabel('zscore');
                title('Peri-Event Time Histogram with Firing Rate');
                xlim(psth_window);
                ylim([-0.5,1.5])
                %             sgtitle([ animal '-day-' num2str(curr_file) ' ' spike_name ])


            end

            saveas(gcf, fullfile(Path,newfolderName,[ animal  '_day_' num2str(curr_file) spike_name '_PSTH.jpg']),'jpg')
            close all
        end
    
%      for curr_cell=1:length(spike_whole)
%             spike_times=spike_whole{curr_cell};
%             spike_name=neuron_files(curr_cell).name(1:end-4);
% 
%   
%    idx=spike_times>event_times{1}'&spike_times<event_times{3}';
% 
%    for curr_trial=1:length(event_times{1})
%        spike_n=spike_times(idx(:,curr_trial))
%        norm_spike{curr_trial}=(spike_n-spike_n(1))/(spike_n(end)-spike_n(1));
%    end


    
    
    end
end
