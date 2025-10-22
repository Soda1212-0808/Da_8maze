clear all
clc

Path=ds.locations.server_data_path;

load_parts=struct;
load_parts.tetrode=true;
load_parts.video_track=true;

animal='DCA3-9';
rec_day='2021-06-13';

ds.load_recoding


time_intervals=cellfun(@(tb,tc) ...
    arrayfun(@(i) position_timelite >= tb(i) & position_timelite <= tc(i), (1:size(tb,1))', ...
    'UniformOutput', false),arm_times{1},arm_times{3},'UniformOutput',false);

%% 计算占用直方图

bin_size=10;
frame_rate=30;
inIntervals=cell2mat(cat(1, time_intervals{:})');
[occupancy_time,x_edges,y_edges]=ds.occupancy(position_re_X,position_re_Y,'timeinterval',inIntervals,'bin_size',10,'fps',30);





X_by_trial=arrayfun(@(trial) position_re_X(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
Y_by_trial=arrayfun(@(trial) position_re_Y(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
X_by_trial_filled=cellfun(@(x)  interp1(find(~isnan(x)),x(~isnan(x)),(1:length(x))',"linear"),X_by_trial, 'UniformOutput', false)
Y_by_trial_filled=cellfun(@(x)  interp1(find(~isnan(x)),x(~isnan(x)),(1:length(x))',"linear"),Y_by_trial, 'UniformOutput', false)
X_resort=cell2mat(cellfun(@(x)  [x;nan],X_by_trial_filled,'UniformOutput',false)');
Y_resort=cell2mat(cellfun(@(x)  [x;nan],Y_by_trial_filled,'UniformOutput',false)');


position_time_by_trial=arrayfun(@(trial) position_timelite(inIntervals(:,trial)) ,1:size(inIntervals,2), 'UniformOutput', false)
position_time_resort=cell2mat(cellfun(@(x)  [x;x(end)+1/frame_rate],position_time_by_trial,'UniformOutput',false)');
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


