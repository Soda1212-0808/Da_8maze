clear all
% 定义包含CSV文件的文件夹路径
Path = 'G:\CA3_rawdata\CA3_2p\data';    % 设置数据存放的文件夹路径
% animals={'1464'};
animals={'1646','1306','1307','1309','1311','1312','1974','1976'};

curr_animal=2
animal=animals{curr_animal};

newfolderName = 'buffer_image';
if exist(fullfile(Path,animal,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,animal,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end

contents = dir(fullfile(Path ,animal));
% 获取所有子文件夹的名称
recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

bufferfolderName = 'bufferFile';
if exist(fullfile(Path, animal ,bufferfolderName), 'dir') ~= 7
    mkdir(fullfile(Path, animal , bufferfolderName));
    disp(['Folder "', bufferfolderName, '" created.']);
end

all_data_path=load(fullfile(Path,animal,'merged file','merged_mice_path.mat'));


X_position=cellfun(@(x)table2array(x(:,2))+100,all_data_path.animal_path  ,'UniformOutput',false);
Y_position=cellfun(@(x)table2array(x(:,3))+100,all_data_path.animal_path  ,'UniformOutput',false);

%% 前处理迷宫轨迹采集错误的数据

%获取MP4文件
mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.mp4'));
if isempty(mp4Files)
    mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.avi'));
end
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
    for k = 6:numPolygons
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
all_data_X_filter_speed=cell(length(X_position),1);
all_data_Y_filter_speed=cell(length(X_position),1);
occupancy_time_all=cell(length(X_position),1);
% x_edges=cell(length(X_position),1);
% y_edges=cell(length(X_position),1);
 x_edges = 0:bin_size:1200;
 y_edges = 0:bin_size:1200;

all_data_speed=cell(length(X_position),1);
figure('Name',animal,'Position',[50 50 800 800]);
%设置速度阈值
speed_threshold=0;
for curr_day=1:length(X_position)
    X_interp=X_position{curr_day};
    Y_interp=Y_position{curr_day};
    nexttile
    plot(X_interp,Y_interp)
    title(['raw day ' num2str(curr_day) ])

    % 找到 mask 中对应位置为 0 的索引
    mask=BW1{1};
    invalid_indices = mask(sub2ind(size(mask), fix(Y_interp+1), fix(X_interp+1))) == 0;
    % 将对应的 X_interp 和 Y_interp 位置设为 NaN
    X_interp(invalid_indices) = NaN;
    Y_interp(invalid_indices) = NaN;
    nexttile
        plot(X_interp,Y_interp)

    window_size = 5;
    distances_X = abs(X_interp - movmean(X_interp, [window_size window_size], 'omitnan'));
    distances_Y = abs(Y_interp - movmean(Y_interp, [window_size window_size], 'omitnan'));
    % 计算欧几里得距离
    distances = sqrt(distances_X.^2 + distances_Y.^2);
    nexttile;
    plot(distances)
    threshold=40;
    errors = distances > threshold;
    X_interp(errors)=NaN;
    Y_interp(errors)=NaN;

    nexttile
    plot(X_interp,Y_interp);xlim([0 1200]);ylim([0 1200]);
    title(['processed day' num2str(curr_day)])
    speed = [0 ;sqrt(diff(X_interp).^2 + diff(X_interp).^2)];
    nexttile
    plot(speed);
    title('speed')

    X_filter_speed=X_interp; X_filter_speed(speed<speed_threshold)=NaN;
    Y_filter_speed=Y_interp; Y_filter_speed(speed<speed_threshold)=NaN;


    %%网格分辨率
    bin_size=10;
    % data_path.time=(data_path.scorer)/framerate;
    % position_time= data_path.time;
    frame_sampling_rate=30;
    % 定义网格的边界和分辨率
    % x_edges{curr_day} = min(X_filter_speed):bin_size:max(X_filter_speed);
    % y_edges{curr_day} = min(Y_filter_speed):bin_size:max(Y_filter_speed);
   
    % 计算占用直方图
    %         occupancy_map = histcounts2(X_filter_speed, Y_filter_speed, x_edges{curr_day}, y_edges{curr_day});
    % 计算占用时间（假设位置数据的采样率为 position_sampling_rate）
    %         occupancy_time = occupancy_map * (1 / frame_sampling_rate);

    all_data_X_filter_speed{curr_day}=X_filter_speed;
    all_data_Y_filter_speed{curr_day}=Y_filter_speed;
    all_data_speed{curr_day}=speed;
    %         occupancy_time_all{curr_day}=occupancy_time;
end



%% 计算perievent
all_data_event=load(fullfile(Path,animal,'merged file','merged_mice_behavior_timepoint.mat'));
all_data_match=load(fullfile(Path,animal,'merged file','merged_mice_cell_timepoint.mat'));

all_data_event.all_event_frame = cellfun(@(A,B) [round(interp1(table2array(A(:, 3)), (1:size(A,1))',B(:,1:5), 'linear', 'extrap')) B(:,6:7)],...
    all_data_match.animal_timepoint, all_data_event.all_event_timepoint','UniformOutput',false);


track_sample_arm=BW1{2};
track_reward_arm= BW1{3}+BW1{4}+BW1{5}+BW1{6}+BW1{7}+BW1{8}+BW1{9};

% figure;imagesc(track_reward_arm)
% 提取上述四个event的时间点

sample_begin =cell(length(all_data_event.all_event_frame),1);
choice_begin =cell(length(all_data_event.all_event_frame),1);
sample_arm =cell(length(all_data_event.all_event_frame),1);
choice_arm =cell(length(all_data_event.all_event_frame),1);

for curr_day=1:length(all_data_event.all_event_frame)

    sample_begin{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);
    choice_begin{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);

    sample_arm{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);
    choice_arm{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);


    figure('Position',[5 5 1800 1000]);

    for curr_trial=1:size(all_data_event.all_event_frame{curr_day},1)

        % sample begin time point
        sample_range=all_data_event.all_event_frame{curr_day}(curr_trial,[1,3]);
        sample_x1=all_data_X_filter_speed{curr_day}(sample_range(1):sample_range(2));
        sample_y1=all_data_Y_filter_speed{curr_day}(sample_range(1):sample_range(2));
        sample_in_target_area_1=nan(length(sample_x1),1);

        sample_in_target_area_1(~isnan(sample_x1)) = track_sample_arm(sub2ind(size(track_sample_arm), round(sample_y1(~isnan(sample_x1))), round(sample_x1(~isnan(sample_x1)))));

        if ~isempty(find(sample_in_target_area_1==1,1,'last'))
            sample_begin{curr_day}(curr_trial)=sample_range(1)-1+ find(sample_in_target_area_1==1,1,'last');
        else  sample_begin{curr_day}(curr_trial)=nan;
        end

        sample_in_target_area_2=nan(length(sample_x1),1);
        sample_in_target_area_2(~isnan(sample_x1)) = track_reward_arm(sub2ind(size(track_reward_arm), round(sample_y1(~isnan(sample_x1))), round(sample_x1(~isnan(sample_x1)))));

        if ~isempty(find(sample_in_target_area_2==1,1,'first'))
            sample_arm{curr_day}(curr_trial)=sample_range(1)-1+ find(sample_in_target_area_2==1,1,'first');
        else  sample_arm{curr_day}(curr_trial)=nan;
        end


        nexttile
        hold on
        imagesc(track_sample_arm+track_reward_arm);
        xlim([0 1200]);  ylim([0 1200]);
        axis image off;
        colormap( ap.colormap('WK'));clim([0 10])
        plot(round(sample_x1), round(sample_y1),'LineWidth',2,'Color','r')
        scatter(sample_x1(find(sample_in_target_area_1==1,1,'last')),sample_y1(find(sample_in_target_area_1==1,1,'last')),20,'blue','filled')
        scatter(sample_x1(find(sample_in_target_area_2==1,1,'first')),sample_y1(find(sample_in_target_area_2==1,1,'first')),20,'yellow','filled')
        % sample begin time point
        title(['sample ' num2str(curr_trial)])




        choice_range=[all_data_event.all_event_frame{curr_day}(curr_trial,4) all_data_event.all_event_frame{curr_day}(curr_trial,5)+100];
        choice_x1=all_data_X_filter_speed{curr_day}(choice_range(1):choice_range(2));
        choice_y1=all_data_Y_filter_speed{curr_day}(choice_range(1):choice_range(2));
        choice_in_target_area_1=nan(length(choice_x1),1);
        choice_in_target_area_1(~isnan(choice_x1)) = track_sample_arm(sub2ind(size(track_sample_arm), round(choice_y1(~isnan(choice_x1))), round(choice_x1(~isnan(choice_x1)))));

        if ~isempty(find(choice_in_target_area_1==1,1,'last'))
            choice_begin{curr_day}(curr_trial)=choice_range(1)-1+ find(choice_in_target_area_1==1,1,'last');
        else  choice_begin{curr_day}(curr_trial)=nan;
        end

        choice_in_target_area_2=nan(length(choice_x1),1);
        choice_in_target_area_2(~isnan(choice_x1)) = track_reward_arm(sub2ind(size(track_reward_arm), round(choice_y1(~isnan(choice_x1))), round(choice_x1(~isnan(choice_x1)))));
        if ~isempty(find(choice_in_target_area_2==1,1,'first'))
            choice_arm{curr_day}(curr_trial)=choice_range(1)-1+ find(choice_in_target_area_2==1,1,'first');
        else  choice_arm{curr_day}(curr_trial)=nan;
        end

        nexttile
        hold on
        imagesc(track_sample_arm+track_reward_arm)
        xlim([0 1200]);  ylim([0 1200]);
        axis image off;
        colormap( ap.colormap('WK'));clim([0 10])

        plot(round(choice_x1), round(choice_y1),'LineWidth',2,'Color','r')
        scatter(choice_x1(find(choice_in_target_area_1==1,1,'last')),choice_y1(find(choice_in_target_area_1==1,1,'last')),20,"green",'filled')
        scatter(choice_x1(find(choice_in_target_area_2==1,1,'first')),choice_y1(find(choice_in_target_area_2==1,1,'first')),20,'yellow','filled')
        title(['choice ' num2str(curr_trial)])

        drawnow
        % sample begin time point

    end
    sgtitle(['day ' num2str(curr_day-1)])
    saveas(gcf,fullfile(Path, animal ,bufferfolderName,['signle trial visuaization ' num2str(curr_day-1) '.jpg']),'jpeg')


end

newfile_name=fullfile(Path,animal,'bufferfile', [animal 'event_data.xlsx']);

% 保存到 Excel，不同工作表名称对应变量名

for i = 1:5
    %     writematrix(data{i}, filename, 'Sheet', sheetName); % 写入 Excel
    buffer_data=[sample_begin{i} sample_arm{i} choice_begin{i} choice_arm{i}]
    writematrix(buffer_data, newfile_name, 'Sheet', ['day' num2str(i-1)]); % B 为工作表名称

end

