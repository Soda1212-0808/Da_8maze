clear all
% 定义包含CSV文件的文件夹路径


Path = 'D:\SDdata\data_2p_1';    % 设置数据存放的文件夹路径
animals={'1464'};
animal='1464';
contents = dir(fullfile(Path ,animal));
% 获取所有子文件夹的名称
recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};


bufferfolderName = 'bufferFile';
if exist(fullfile(Path, animal ,bufferfolderName), 'dir') ~= 7
    mkdir(fullfile(Path, animal , bufferfolderName));
    disp(['Folder "', bufferfolderName, '" created.']);
end


all_data_path=load(fullfile(Path,animal,'merged file','merged_mice_path.mat'));
data_path=all_data_path.cell_animal_path{1};

%% load Ca2+数据
data_imaging=load(fullfile(Path,animal,'data_2p_cell\cell_video_align\suite2p\plane0','Fall.mat'));
data_imgaing_curr=data_imaging.spks(:,1:size(data_path,1));


%% 前处理迷宫轨迹采集错误的数据
X=table2array(data_path(:,2))+100;
Y=table2array(data_path(:,9))+100;

%获取MP4文件
mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.mp4'));
mp4FilePaths=fullfile({mp4Files.folder}, {mp4Files.name});
% 创建一个 VideoReader 对象
v = VideoReader(mp4FilePaths{1});
framerate=v.framerate;

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
if exist(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'))~=2
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
saveas(gcf,fullfile(Path, animal ,bufferfolderName,'grab_picture.jpg'),'jpeg')
save(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'),'BW1','-mat')
% save(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'),'BW1','recordedFrameCount','-mat')

else
load(fullfile(Path, animal ,bufferfolderName,'grab_picture.mat'))
end


X_interp=X;
Y_interp=Y;
% 找到 mask 中对应位置为 0 的索引
mask=BW1{1};
invalid_indices = mask(sub2ind(size(mask), fix(Y_interp+1), fix(X_interp+1))) == 0;
% 将对应的 X_interp 和 Y_interp 位置设为 NaN
X_interp(invalid_indices) = NaN;
Y_interp(invalid_indices) = NaN;

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


%%网格分辨率
bin_size=10;
% data_path.time=(data_path.scorer)/framerate;
% position_time= data_path.time;
frame_sampling_rate=30;
% 定义网格的边界和分辨率
x_edges = min(X_filter_speed):bin_size:max(X_filter_speed);
y_edges = min(Y_filter_speed):bin_size:max(Y_filter_speed);

% 计算占用直方图
occupancy_map = histcounts2(X_filter_speed, Y_filter_speed, x_edges, y_edges);
% 计算占用时间（假设位置数据的采样率为 position_sampling_rate）
occupancy_time = occupancy_map * (1 / frame_sampling_rate);

figure('Position',[50 50 1600 800]);
colormap('jet')
for curr_cell=1:size(data_imgaing_curr,1)
    spike_times=spike_whole{curr_cell};
    % 获取每个发放事件对应的位置索引
  
    % 计算发放直方图
    spike_map = histcounts2(X_filter_speed, Y_filter_speed, x_edges, y_edges);
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
    smoothed_rate_map=flipud(smoothed_rate_map);
    nexttile(curr_cell);
    imagesc(x_edges, y_edges, smoothed_rate_map); axis image off;
    clim([0 nanmax(smoothed_rate_map(:))])
    colorbar
    formatted_value = sprintf('%.1f', round(spike_freq(curr_cell),1));
    modified_string = strrep(neuron_files(curr_cell).name(1:end-4), '_', '-');

    title([modified_string ': ' formatted_value 'Hz'])
end
sgtitle([ animal '-day-' num2str(curr_file) ])
saveas(gcf, fullfile(Path,[ animal '_day_' num2str(curr_file) 'rate_map.jpg']),'jpg')
% close all



