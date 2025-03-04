Path='G:\CA3_rawdata\CA3_2p\data\1646\data_path_DLC\';
% 读取视频文件
all_mp4Files = dir(fullfile(Path, '*labeled.mp4'))';
    outputFile=fullfile(all_mp4Files(1).folder,'new', 'dd.mp4');
videoWriter = VideoWriter(outputFile, 'MPEG-4');
open(videoWriter);
frameIndex = 1;
for curr_file =1:length(all_mp4Files)
    inputFile = fullfile(all_mp4Files(curr_file).folder, all_mp4Files(curr_file).name);

% 创建 VideoReader 和 VideoWriter 对象
videoReader = VideoReader(inputFile);

% 设置文本属性
fontSize = 30; % 字体大小
textColor = [0, 255, 0]; % 文字颜色（红色）
position = [10, 10]; % 左上角的位置 (x, y)

% 逐帧处理

while hasFrame(videoReader)
    % 读取当前帧
    frame = readFrame(videoReader);
    
    % 将帧号写入帧图像
    frame = insertText(frame, position, sprintf('Frame: %d', frameIndex), ...
        'FontSize', fontSize, 'TextColor', textColor, 'BoxOpacity', 0);
    
    % 写入到输出视频
    writeVideo(videoWriter, frame);
    
    % 更新帧索引
    frameIndex = frameIndex + 1;
end


end
% 关闭 VideoWriter
close(videoWriter);

disp('处理完成！输出文件为: ' + string(outputFile));