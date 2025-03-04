% 读取视频
videoReader = VideoReader('G:\CA3_rawdata\CA3_2p\B\output_video.mp4');
videoWriter = VideoWriter('G:\CA3_rawdata\CA3_2p\B\output_video1.mp4', 'MPEG-4');
open(videoWriter);

while hasFrame(videoReader)
    frame = readFrame(videoReader); % 读取一帧

    % 转换为灰度（可选）
    % frame = rgb2gray(frame);

    % 调整对比度和亮度
    adjustedFrame1 = imadjust(frame, stretchlim(frame, [0.5 0.95]), [0 1]);
     adjustedFrame = rescale(adjustedFrame1, 0.5, 0.99); % 调整到新的亮度范围

    % 写入新视频
    writeVideo(videoWriter, adjustedFrame);
end

close(videoWriter);
disp('视频处理完成');
