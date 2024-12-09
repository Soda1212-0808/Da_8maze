clear all

Path = 'H:\CA3_reprocess\each mice';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};

% Path = 'H:\CA1_8Maze';    % 设置数据存放的文件夹路径
% animals={'CA1-9'};

for curr_animal=1:length(animals)
% 获取文件夹中的所有内容
% animal='CA1-9';
animal=animals{curr_animal};
contents = dir(fullfile(Path ,animal));
% 获取所有子文件夹的名称
recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

%选择文件
curr_file=1;
% curr_file=1:length(recording_files)

%%
bufferfolderName = 'bufferFile';
if exist(fullfile(Path, animal , recording_files{curr_file},bufferfolderName)) ~= 7
    mkdir(fullfile(Path, animal , recording_files{curr_file},bufferfolderName));
    disp(['Folder "', bufferfolderName, '" created.']);
end


%获取MP4文件
mp4Files=dir(fullfile(Path, animal , recording_files{curr_file} ,'*.mp4'));
mp4FilePaths=fullfile({mp4Files.folder}, {mp4Files.name});
% 创建一个 VideoReader 对象
v = VideoReader(mp4FilePaths{1});
framerate=v.framerate;

neuron_files=dir(fullfile(Path, animal , recording_files{curr_file} ,'*.t64'));
wv_files=dir(fullfile(Path, animal , recording_files{curr_file} ,'*wv.mat'));

% curr_csv=dir(fullfile(Path,animal,recording_files{curr_file},'*filtered.csv'));
curr_csv=dir(fullfile(Path,animal,recording_files{curr_file},'*.csv'));

data_path=readtable([Path '\' animal '\' recording_files{curr_file} '\' recording_files{curr_file} '-' animal 'DLC_resnet50_old-mazeJul21shuffle1_100000_filtered.csv']);
data_event=readtable([Path '\' animal '\' recording_files{curr_file} '\' recording_files{curr_file} '-' animal '-01Dat_data_m1.csv']);

data_path.time=(data_path.scorer+1)/framerate;

spike_whole = arrayfun(@(f) readmclusttfile(fullfile(f.folder, f.name))'/10000, ...
    neuron_files, 'UniformOutput', false);
spike_freq=cell2mat(cellfun(@(x) length(x)/x(end),spike_whole,'UniformOutput',false));

% plot(spike_wv{1, 1}.xrange ,spike_wv{1, 1}.mWV  )
%% 计算峰谷宽
spike_wv = arrayfun(@(f) load(fullfile(f.folder, f.name)), ...
    wv_files, 'UniformOutput', false);

% 定义匿名函数，计算每列的最大值和最小值的差值
calcDifferences = @(mWV) max(mWV) - min(mWV);

% 使用 cellfun 对每个 cell 应用匿名函数，计算每列的差值
differences = cellfun(@(x) calcDifferences(x.mWV), spike_wv, 'UniformOutput', false);

% 找到每个 cell 中差值最大的那一列
[~, maxDiffCols] = max(reshape(cell2mat(differences),length(differences),4), [], 2);

% 定义匿名函数，计算每个 cell 中最大值与最小值之间的距离
calcDistance = @(x, col) abs(find(x.mWV(:, col) == max(x.mWV(:, col)), 1) - ...
    find(x.mWV(:, col) == min(x.mWV(:, col)), 1));
% 使用 cellfun 计算每个 cell 中最大值与最小值之间的距离
wv_distances = arrayfun(@(i) calcDistance(spike_wv{i}, maxDiffCols(i)), 1:numel(spike_wv));

figure;
scatter(spike_freq,wv_distances)
xlabel('firing rate')
ylabel('peak valley distance')


%% 前处理迷宫轨迹采集错误的数据
X=table2array(data_path(:,2))+100;
Y=table2array(data_path(:,9))+100;
 plot(X,Y)

% v.NumFrames
% 读取第一帧
firstFrame = readFrame(v);
[rows, cols, channels] = size(firstFrame);
% 创建一个新图像，尺寸比原始图像大200像素（上下左右各100像素）
newRows = rows + 200;
newCols = cols + 200;
paddedFrame = uint8(255 * ones(newRows, newCols, channels)); % 白色背景
% 将原始图像复制到新图像的中心
paddedFrame(101:100+rows, 101:100+cols, :) = firstFrame;

%% 若之前未绘制目标区域并保存文件，在轨迹图像上绘制目标区域
recordedFrameCount=1;

if ~exist(fullfile(Path, animal , recording_files{curr_file},bufferfolderName,'grab_picture.jpg'), 'file')
display_next_frame_on_scroll(mp4FilePaths{1})

figure;
% 显示填充后的第一帧
imshow(paddedFrame);
hold on
scatter(X,Y)
% plot(X,Y)
% mask=roipoly;
% 指定需要绘制的多边形区域数量
numPolygons = 7; % 你可以根据需要改变此值
BW1=cell(1,numPolygons );
for k = 1:numPolygons
    % 绘制多边形区域
    BW = roipoly;
    BW1{k}=BW;
    % 将当前多边形区域添加到组合掩码中
    %     BW_combined = BW_combined | BW;

    % 显示当前多边形区域的边界
    boundary = bwboundaries(BW);
    plot(boundary{1}(:,2), boundary{1}(:,1), 'LineWidth', 2);
    mark_x = boundary{1}(1, 2);
    mark_y = boundary{1}(1, 1);
    text(mark_x, mark_y, num2str(k-1), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');

end
saveas(gcf,fullfile(Path, animal , recording_files{curr_file} ,bufferfolderName,'grab_picture.jpg'),'jpeg')

save(fullfile(Path, animal , recording_files{curr_file} ,bufferfolderName,'grab_picture.mat'),'BW1','recordedFrameCount','-mat')
% save(fullfile(Path, animal , recording_files{curr_file} ,bufferfolderName,'grab_picture.mat'),'BW1','-mat')

else
load(fullfile(Path, animal , recording_files{curr_file} ,bufferfolderName,'grab_picture.mat'))
end

X_interp=X(recordedFrameCount:end);
Y_interp=Y(recordedFrameCount:end);

% 找到 mask 中对应位置为 0 的索引
mask=BW1{1};
invalid_indices = mask(sub2ind(size(mask), fix(Y_interp+1), fix(X_interp+1))) == 0;
% 将对应的 X_interp 和 Y_interp 位置设为 NaN
X_interp(invalid_indices) = NaN;
Y_interp(invalid_indices) = NaN;


%%
window_size = 5;

distances_X = abs(X_interp - movmean(X_interp, [window_size window_size], 'omitnan'));
distances_Y = abs(Y_interp - movmean(Y_interp, [window_size window_size], 'omitnan'));

% 计算欧几里得距离
distances = sqrt(distances_X.^2 + distances_Y.^2);

figure;
plot(distances)
threshold=40;
errors = distances > threshold;
X_interp(errors)=NaN;
Y_interp(errors)=NaN;

figure;
plot(X_interp,Y_interp);

speed = [0 ;sqrt(diff(X_interp).^2 + diff(X_interp).^2)];

X_filter_speed=X_interp; X_filter_speed(speed<0.5)=NaN;
Y_filter_speed=Y_interp; Y_filter_speed(speed<0.5)=NaN;


% % 定义高斯滤波器的标准差（窗口大小）
% sigma = 2; % 根据你的数据和需要进行调整
% % 对位置数据进行高斯平滑
% X_interp = imgaussfilt(X_interp, sigma);
% Y_interp = imgaussfilt(Y_interp, sigma);

%%网格分辨率
bin_size=10;
data_path.time=(data_path.scorer+1-recordedFrameCount)/framerate;

position_time= data_path.time(recordedFrameCount:end);

inIntervals = any(position_time >= (table2array(data_event( :,4))-0)' & position_time <= table2array(data_event(:,5))', 2);
% inIntervals=(position_time>0);

position_time2=position_time(inIntervals);
X_filter_speed2=X_filter_speed(inIntervals);
Y_filter_speed2=Y_filter_speed(inIntervals);

frame_sampling_rate=30;
% 定义网格的边界和分辨率
x_edges = min(X_filter_speed2):bin_size:max(X_filter_speed2);
y_edges = min(Y_filter_speed2):bin_size:max(Y_filter_speed2);

% 计算占用直方图
occupancy_map = histcounts2(X_filter_speed2, Y_filter_speed2, x_edges, y_edges);
% 计算占用时间（假设位置数据的采样率为 position_sampling_rate）
occupancy_time = occupancy_map * (1 / frame_sampling_rate);

figure('Position',[50 50 1600 800]);
colormap('jet')
for curr_cell=1:length(spike_whole)
    spike_times=spike_whole{curr_cell};
    % 获取每个发放事件对应的位置索引
    spike_x = interp1(position_time2, X_filter_speed2, spike_times);
    spike_y = interp1(position_time2, Y_filter_speed2, spike_times);

    % 计算发放直方图
    spike_map = histcounts2(spike_x, spike_y, x_edges, y_edges);
    % 计算发放速率地图
    rate_map = spike_map ./ occupancy_time;
    rate_map(isnan(rate_map)) = 0; % 将NaN值（由于0占用时间导致的）设为0
    rate_map(isinf(rate_map)) = 0; % 将NaN值（由于0占用时间导致的）设为0

    smooth_sigma = 2; % 根据需要调整
    % 对发放速率地图进行高斯平滑
    smoothed_rate_map = imgaussfilt(rate_map, smooth_sigma);

    % 计算平均发放速率

    mean_rate = mean(rate_map(rate_map>0));
    % 识别发放场
    threshold = 1.5 * mean_rate; % 设置阈值为平均发放速率的两倍
    firing_field = rate_map > threshold;
%     smoothed_rate_map=flipud(smoothed_rate_map);
    nexttile(curr_cell);
   
    imagesc(x_edges, y_edges, smoothed_rate_map); axis image off;
   
    clim([0 nanmax(smoothed_rate_map(:))])
%     colorbar
    formatted_value = sprintf('%.1f', round(spike_freq(curr_cell),1));
    modified_string = strrep(neuron_files(curr_cell).name(1:end-4), '_', '-');

    title([modified_string ': ' formatted_value 'Hz'])
end
sgtitle([ animal '-day-' num2str(curr_file) ])

saveas(gcf, fullfile(Path,[ animal '_day_' num2str(curr_file) 'rate_map.jpg']),'jpg')
% close all





end

