clear all
% 定义包含CSV文件的文件夹路径
Path = 'G:\CA3_rawdata\CA3_2p\data';    % 设置数据存放的文件夹路径
% animals={'1464'};
animals={'1306','1307','1309','1311','1312','1646','1974','1976'};

for curr_animal=6
    animal=animals{curr_animal};

newfolderName = 'buffer_image';
if exist(fullfile(Path,animal,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,animal,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end


% animal='1306';
contents = dir(fullfile(Path ,animal));
% 获取所有子文件夹的名称
recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};


bufferfolderName = 'bufferFile';
if exist(fullfile(Path, animal ,bufferfolderName), 'dir') ~= 7
    mkdir(fullfile(Path, animal , bufferfolderName));
    disp(['Folder "', bufferfolderName, '" created.']);
end


all_data_path=load(fullfile(Path,animal,'merged file','merged_mice_path.mat'));
all_data_match=load(fullfile(Path,animal,'merged file','merged_mice_cell_timepoint.mat'));


file_name=fieldnames(all_data_match.animal_match_table);
lengths = structfun(@(f)  size(f,1)    ,all_data_match.animal_match_table, 'UniformOutput', false);
cut_edge=cell2mat(struct2cell(lengths));

% 仅针对1464
cut_edge(end)=cut_edge(end)-1;


%% load Ca2+数据
data_imaging=load(fullfile(Path,animal,'data_2p_cell\cell_video_align\suite2p\plane0','Fall.mat'));
startIndices = [1; cumsum(cut_edge(1:end-1)) + 1]; % 起始列索引
endIndices = cumsum(cut_edge);                    % 结束列索引
% 使用数组操作生成子矩阵的 cell 数组
subMatrices_imaging = arrayfun(@(s, e) data_imaging.spks(:, s:e), startIndices, endIndices, 'UniformOutput', false);

data_imgaing_all=cell(5,1);
data_path=cell(5,1);
day_name={'day 0','day 1','day 2','day 3','day 4'};
for curr_day=1:5
% curr_day=4

match_id=all_data_match.animal_match_table.(file_name{curr_day});
data_path{curr_day}=all_data_path.cell_animal_path{curr_day}(cell2mat(match_id(:,3)),:);
data_imgaing_all{curr_day}=subMatrices_imaging{curr_day};
end

%% 前处理迷宫轨迹采集错误的数据

X_position=cellfun(@(x)table2array(x(:,2))+100,data_path,'UniformOutput',false);
Y_position=cellfun(@(x)table2array(x(:,9))+100,data_path,'UniformOutput',false);

%获取MP4文件
mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.mp4'));
%  mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.avi'));
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

% 若之前未绘制目标区域并保存文件，在轨迹图像上绘制目标区域
if exist(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'))~=2
display_next_frame_on_scroll(mp4FilePaths{1})
figure;
% 显示填充后的第一帧
imshow(paddedFrame);
hold on
cellfun(@(x,y) scatter(x,y),X_position,Y_position,'UniformOutput',false)

% plot(X,Y)
% mask=roipoly;
% 指定需要绘制的多边形区域数量
numPolygons = 9; % 你可以根据需要改变此值
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
%% 提取轨迹
X_filter_speed_all=cell(length(X_position),1);
Y_filter_speed_all=cell(length(X_position),1);
occupancy_time_all=cell(length(X_position),1);
x_edges=cell(length(X_position),1);
y_edges=cell(length(X_position),1);
for curr_day=1:length(X_position)
X_interp=X_position{curr_day};
Y_interp=Y_position{curr_day};
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
plot(X_interp,Y_interp);xlim([0 1200]);ylim([0 1200]);
speed = [0 ;sqrt(diff(X_interp).^2 + diff(X_interp).^2)];
X_filter_speed=X_interp; X_filter_speed(speed<0.5)=NaN;
Y_filter_speed=Y_interp; Y_filter_speed(speed<0.5)=NaN;


%%网格分辨率
bin_size=10;
% data_path.time=(data_path.scorer)/framerate;
% position_time= data_path.time;
frame_sampling_rate=30;
% 定义网格的边界和分辨率
% x_edges{curr_day} = min(X_filter_speed):bin_size:max(X_filter_speed);
% y_edges{curr_day} = min(Y_filter_speed):bin_size:max(Y_filter_speed);
x_edges{curr_day} = 0:bin_size:1200;
y_edges{curr_day} = 0:bin_size:1200;

% 计算占用直方图
occupancy_map = histcounts2(X_filter_speed, Y_filter_speed, x_edges{curr_day}, y_edges{curr_day});
% 计算占用时间（假设位置数据的采样率为 position_sampling_rate）
occupancy_time = occupancy_map * (1 / frame_sampling_rate);

X_filter_speed_all{curr_day}=X_filter_speed;
Y_filter_speed_all{curr_day}=Y_filter_speed;
occupancy_time_all{curr_day}=occupancy_time;
end


%%  计算 神经元在位置上的放电热区图   firing map
for curr_cell=1:size(data_imgaing_all{1},1)

    figure('Position',[50 50 800 200]);
    colormap('jet')
    for curr_day=1:length(data_imgaing_all)
        spike_times=data_imgaing_all{curr_day}(curr_cell,:)';
        % 获取每个发放事件对应的位置索引
        spike_x{curr_day} = X_filter_speed_all{curr_day}(spike_times>50);
        spike_y{curr_day} = Y_filter_speed_all{curr_day}(spike_times>50);

        % 计算发放直方图
        spike_map{curr_day} = histcounts2(spike_x{curr_day}, spike_y{curr_day}, x_edges{curr_day}, y_edges{curr_day});
        % 计算发放速率地图
        rate_map{curr_day} = spike_map{curr_day} ./ occupancy_time_all{curr_day};
        rate_map{curr_day}(isnan(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0
        rate_map{curr_day}(isinf(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0

        smooth_sigma = 2; % 根据需要调整
        % 对发放速率地图进行高斯平滑
        rate_map_smoothed{curr_day} = imgaussfilt(rate_map{curr_day}, smooth_sigma);

        % 计算平均发放速率

        mean_rate{curr_day} = mean(rate_map{curr_day}(rate_map{curr_day}>0));
        % 识别发放场
        threshold = 1.5 * mean_rate{curr_day}; % 设置阈值为平均发放速率的两倍
        firing_field{curr_day} = rate_map{curr_day} > threshold;
        %%是否反转图像
        %     smoothed_rate_map=flipud(smoothed_rate_map);
        rate_map_smoothed{curr_day}=rate_map_smoothed{curr_day};

        nexttile(curr_day);

        imagesc(x_edges{curr_day}, y_edges{curr_day}, rate_map_smoothed{curr_day});
        axis image off;
        clim([0 max(rate_map_smoothed{1}(:))])
        xlim([0 1200]);     ylim([0 1200]);

        %     colorbar
        %     formatted_value = sprintf('%.1f', round(spike_freq(curr_cell),1));
        %     modified_string = strrep(neuron_files(curr_cell).name(1:end-4), '_', '-');
        %
        %     title([modified_string ': ' formatted_value 'Hz'])
        title(day_name{curr_day})
        drawnow
    end

    colorbar
    sgtitle([ animal '-neuron-' num2str(curr_cell) ])

saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal 'neuron' num2str(curr_cell) 'rate_map.jpg']),'jpg')

end
% close all

%% 计算perievent
all_data_event=load(fullfile(Path,animal,'merged file','merged_mice_behavior_timepoint.mat'));

all_data_event.all_event_frame = cellfun(@(A,B) round(interp1(table2array(A(:, 4)), (1:size(A,1))',B, 'linear', 'extrap')),...
    all_data_match.cell_timepoint_cell, all_data_event.all_event_timepoint','UniformOutput',false);




end
