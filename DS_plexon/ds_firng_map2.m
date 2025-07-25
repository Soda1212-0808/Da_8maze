clear all

% Path = 'H:\Telc-CA3\each mice';
% animals={'DCA3-22(LCA3-Telc)','DCA3-23(LCA3-Telc)','DCA3-24(RCA3-Telc)','DCA3-25(LCA3-Telc)',...
%     'DCA3-26(LCA3-Telc)','DCA3-28(RCA3-Telc)','DCA3-29(RCA3-Telc)'};

Path = 'E:\SDdata\WT\each mice';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};

for curr_animal=7:length(animals)
animal=animals{curr_animal};
contents = dir(fullfile(Path ,animal));
recording_days = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

%选择文件
curr_day=1;
rec_day=recording_days{curr_day};

% load spikes
ds.load_spikes

%
ds.load_video_track
ds.load_events


bin_size=10;
position_time= data_path.time(recordedFrameCount:end);

% sample period or choice period
inIntervals_sample = position_time >= (data_event( :,3)-1)' & position_time <= (data_event(:,5)+1)';
inIntervals_choice = position_time >= (data_event( :,7)-1)' & position_time <= (data_event(:,9)+1)';
inIntervals=[inIntervals_sample inIntervals_choice];

position_time_by_trial=arrayfun(@(trial) position_time(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
X_by_trial=arrayfun(@(trial) X_filter(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
Y_by_trial=arrayfun(@(trial) Y_filter(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
X_by_trial_filled=cellfun(@(x)  interp1(find(~isnan(x)),x(~isnan(x)),(1:length(x))',"linear"),X_by_trial, 'UniformOutput', false)
Y_by_trial_filled=cellfun(@(x)  interp1(find(~isnan(x)),x(~isnan(x)),(1:length(x))',"linear"),Y_by_trial, 'UniformOutput', false)

% figure;
% for curr_trial=1:length(Y_by_trial)
% nexttile
% imagesc(BW1{1})
% hold on
% 
%  plot(X_by_trial_filled{curr_trial},Y_by_trial_filled{curr_trial},'K')
% 
%  plot(X_by_trial{curr_trial},Y_by_trial{curr_trial},'r')
%  xlim([100 700])
%  ylim([0 600])
%  axis square
% end

position_time_resort=cell2mat(cellfun(@(x)  [x;x(end)+1/framerate],position_time_by_trial,'UniformOutput',false)');
X_resort=cell2mat(cellfun(@(x)  [x;nan],X_by_trial_filled,'UniformOutput',false)');
Y_resort=cell2mat(cellfun(@(x)  [x;nan],Y_by_trial_filled,'UniformOutput',false)');
% figure;
% plot(X_resort,Y_resort)

frame_rate=30;
% 定义网格的边界和分辨率
x_edges = min(X_resort):bin_size:max(X_resort);
y_edges = min(Y_resort):bin_size:max(Y_resort);

% 计算占用直方图
occupancy_map = histcounts2(X_resort, Y_resort, x_edges, y_edges);
occupancy_time = occupancy_map * (1 / frame_rate);


% figure;
% nexttile
% imagesc(occupancy_time)


% figure('Position',[50 50 1600 800]);
% colormap('jet')
for curr_cell=1:length(spikes_all)
    spike_times=spikes_all{curr_cell};
    % 获取每个发放事件对应的位置索引

    spike_x = interp1(position_time_resort, X_resort, spike_times);
    spike_y = interp1(position_time_resort, Y_resort, spike_times);


    % 计算发放直方图
    spike_map = histcounts2(spike_x, spike_y, x_edges, y_edges);
    % 计算发放速率地图
    rate_map = spike_map ./ occupancy_time;


    nan_mask = isnan(rate_map);
    rate_map_nan0 = rate_map;
    rate_map_nan0(nan_mask) = 0; % 将NaN值（由于0占用时间导致的）设为0
    rate_map_nan0(isinf(rate_map_nan0))=0;
    smooth_sigma = 1; % 根据需要调整
    % 对发放速率地图进行高斯平滑
    smoothed_rate_map = imgaussfilt(rate_map_nan0, smooth_sigma);

    %          % identify firing field
    %         mean_rate = nanmean(rate_map_nan0(rate_map_nan0>0));
    %         threshold = 1.5 * mean_rate; % 设置阈值为平均发放速率的两倍
    %         firing_field = rate_map > threshold;




    %     smoothed_rate_map=flipud(smoothed_rate_map);
%     nexttile(curr_cell);
    modified_string = strrep(neuron_files(curr_cell).name(1:end-2), '_', '-');

    figure('Name',modified_string)
    % imagesc(x_edges, y_edges, firing_field); axis image off;
    % clim([0 1])
    imagesc(x_edges, y_edges, smoothed_rate_map); axis image off;
    clim([0 nanmax(smoothed_rate_map(:))])
    colormap(ap.colormap('WK'))

    %     colorbar
    formatted_value = sprintf('%.1f', round(spike_freq(curr_cell),1));

    title([modified_string ': ' formatted_value 'Hz'])
    drawnow
end

sgtitle([ animal '-day-' num2str(curr_day) ])

saveas(gcf, fullfile(Path,[ animal '_day_' num2str(curr_day) 'rate_map.jpg']),'jpg')
% close all





end



