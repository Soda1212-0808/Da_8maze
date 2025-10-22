% load video track

% 获取MP4文件
% mp4_file=dir(fullfile(Path, animal , rec_day, 'video_track' ,'*.mp4'));
mp4_file=dir(fullfile(ds.locations.filename('local',animal,rec_day,'video_track'),'*.AVI'));

% 创建一个 VideoReader 对象
v = VideoReader(fullfile(mp4_file.folder, mp4_file.name));
framerate=v.framerate;

% 获取DLC path 文件
path_file=dir(fullfile(Path,animal, rec_day, 'video_track','*.csv'));
data_path=readtable(fullfile(path_file.folder, path_file.name));
data_path.time=(data_path.scorer+1)/framerate;

%% 前处理迷宫轨迹采集错误的数据
poss_nose=table2array(data_path(:,4));
poss_threshold=0.8;
X=table2array(data_path(:,2))+100;
Y=table2array(data_path(:,3))+100;

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

if ~exist(fullfile(Path, animal , rec_day, 'video_track','grab_picture.jpg'), 'file')

figure;
% 显示填充后的第一帧
imshow(paddedFrame);
hold on
scatter(X,Y)

% 指定需要绘制的多边形区域数量
numPolygons = 10; % 你可以根据需要改变此值
BW1=cell(1,numPolygons );
for k = 1:numPolygons
    % 绘制多边形区域
    BW = roipoly;
    BW1{k}=BW;

    % 显示当前多边形区域的边界
    boundary = bwboundaries(BW);
    plot(boundary{1}(:,2), boundary{1}(:,1), 'LineWidth', 2);
    mark_x = boundary{1}(1, 2);
    mark_y = boundary{1}(1, 1);
    text(mark_x, mark_y, num2str(k-1), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');

end

saveas(gcf,fullfile(Path, animal , rec_day, 'video_track','grab_picture.jpg'),'jpeg')
save(fullfile(Path, animal ,rec_day, 'video_track','grab_picture.mat'),'BW1','recordedFrameCount','-mat')

else
load(fullfile(Path, animal , rec_day, 'video_track','grab_picture.mat'))
end

X_filter=X(recordedFrameCount:end);
Y_filter=Y(recordedFrameCount:end);

% 排除识别到的点在迷宫框图之外的
mask=BW1{1};
invalid_indices = mask(sub2ind(size(mask), fix(Y_filter+1), fix(X_filter+1))) == 0;
X_filter(invalid_indices) = NaN;
Y_filter(invalid_indices) = NaN;

% 排除识别置信度小于0.8的
X_filter(poss_nose<0.8) = NaN;
Y_filter(poss_nose<0.8) = NaN;

% 排除两点距离大于40的
window_size = 5;
distances_X = abs(X_filter - movmean(X_filter, [window_size window_size], 'omitnan'));
distances_Y = abs(Y_filter - movmean(Y_filter, [window_size window_size], 'omitnan'));
distances = sqrt(distances_X.^2 + distances_Y.^2);
threshold=40;
X_filter(distances > threshold)=NaN;
Y_filter(distances > threshold)=NaN;


idx_non_nan=~isnan(X_filter);






