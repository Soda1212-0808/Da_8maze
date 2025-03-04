clc; clear; close all;

% 设定主文件夹路径
mainFolder = 'E:\ZZL\视频预处理\all\bycode'; % 修改为你的实际路径
outputFolder='E:\ZZL\视频预处理\all\processed';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); % 如果不存在，则创建文件夹

end


mp4files = dir(fullfile(mainFolder,'*.mp4'));
% 读取视频文件
videoReader = VideoReader(fullfile(mainFolder, mp4files(1).name));

% 读取第一帧
firstFrame = readFrame(videoReader);

% 交互式绘制不规则遮挡区域（生成 mask）
figure, imshow(firstFrame);
title('请手动绘制遮挡区域，然后双击确认');
h = impoly; % 创建可拖动的多边形
wait(h); % 等待用户完成绘制
mask = createMask(h); % 生成二值掩码
close; % 关闭图像窗口


for curr_file=1:length(mp4files)
    fullPath = fullfile(mainFolder, mp4files(curr_file).name);

% 创建视频写入对象
videoWriter = VideoWriter(fullfile(outputFolder, mp4files(curr_file).name) , 'MPEG-4');
videoWriter.FrameRate = videoReader.FrameRate; % 保持原视频帧率
open(videoWriter);

% 重新打开视频文件，以处理所有帧
videoReader = VideoReader(fullPath);

% 逐帧处理视频
while hasFrame(videoReader)
    frame = readFrame(videoReader); % 读取一帧
    
    % 应用 mask，在不规则区域填充白色
    for c = 1:3 % 对 RGB 三通道操作
        frame(:,:,c) = uint8(mask).* frame(:,:,c)  + uint8(~mask) * 255;
    end


    % 转换为灰度
    grayFrame = rgb2gray(frame);
    
    % 曝光
    exposureAdjustment = 30;  % 曝光调整值，可以根据需要增加或减少
    adjustedExposureFrame = grayFrame + exposureAdjustment;  % 增加亮度来提高曝光

    % 对比度
    adjustedContrastFrame = imadjust(adjustedExposureFrame, [0.3 0.9], [0 1]);

    % 将灰度图转换为 RGB 图像，因为 VideoWriter 需要 RGB 输入
    adjustedFrame = repmat(adjustedContrastFrame, [1, 1, 3]);



    % 写入新的视频
    writeVideo(videoWriter, adjustedFrame);
end

% 关闭视频写入对象
close(videoWriter);
disp('视频处理完成，已保存为 output.mp4');
end

