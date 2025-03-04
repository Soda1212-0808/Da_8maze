clc; clear; close all;

% 1️⃣ 设置输入和输出文件夹
input_folder = 'G:\CA3_rawdata\CA3_2p\data\1646\data_path_DLC\';  % 你的 MP4 文件夹路径
output_folder = 'G:\CA3_rawdata\CA3_2p\data\1646\data_path_DLC\new'; % 输出文件夹

if ~exist(output_folder, 'dir')
    mkdir(output_folder); % 如果输出文件夹不存在，则创建
end

% 2️⃣ 获取所有 MP4 文件
video_files = dir(fullfile(input_folder, '*labeled.mp4'));

% 3️⃣ 按日期分组文件
video_groups = containers.Map(); % 创建一个字典存放日期 -> 视频列表

for i = 1:length(video_files)
    filename = video_files(i).name;
    % 提取日期部分（假设格式为 xxxx-xx-xx）
    date_match = regexp(filename, '\d{4}-\d{2}-\d{2}', 'match');
    
    if ~isempty(date_match)
        video_date = date_match{1}; % 获取匹配的日期字符串
        % 将视频按日期存入字典
        if isKey(video_groups, video_date)
            video_groups(video_date) = [video_groups(video_date), {filename}];
        else
            video_groups(video_date) = {filename};
        end
    end
end

% 4️⃣ 逐天合并视频
date_keys = keys(video_groups);

for k = 1:length(date_keys)
    date_str = date_keys{k};
    video_list = video_groups(date_str);
    fprintf('正在合并 %s 的视频...\n', date_str);

    % 读取第一个视频，获取参数
    first_video = VideoReader(fullfile(input_folder, video_list{1}));
    frame_rate = first_video.FrameRate;  % 统一帧率
    video_height = first_video.Height;   % 统一高度
    video_width = first_video.Width;     % 统一宽度

    % 创建输出文件
    output_file = fullfile(output_folder, [date_str, '_merged.mp4']);
    output_video = VideoWriter(output_file, 'MPEG-4');
    output_video.FrameRate = frame_rate;
    open(output_video);

    frame_count = 0; % 记录帧数

    % 逐个处理当天的所有视频
    for i = 1:length(video_list)
        video_path = fullfile(input_folder, video_list{i});
        video_reader = VideoReader(video_path);
        fprintf('   -> 处理 %s\n', video_list{i});
        
        while hasFrame(video_reader)
            frame = readFrame(video_reader);
            frame_count = frame_count + 1;
            
            % 在左上角绘制帧数
            fontSize = 30; % 字体大小
            textColor = [0, 255, 0]; % 文字颜色（红色）
            position = [10, 10]; % 左上角的位置 (x, y)
            frame = insertText(frame, position, sprintf('Frame: %d', frame_count), ...
                'FontSize', fontSize, 'TextColor', textColor, 'BoxOpacity', 0);
%             position = [10, 20]; % 文字位置
%             text_str = sprintf('Frame: %d', frame_count);
%             frame = insertText(frame, position, text_str, 'FontSize', 18, 'BoxColor', 'black', 'TextColor', 'white');
   
            % 统一分辨率
            if size(frame, 1) ~= video_height || size(frame, 2) ~= video_width
                frame = imresize(frame, [video_height, video_width]);
            end
            
            writeVideo(output_video, frame);
        end
    end

    % 关闭视频写入器
    close(output_video);
    fprintf('✅ 合并完成: %s\n', output_file);
end

disp(' 所有视频合并完成！');
