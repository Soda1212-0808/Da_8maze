% 读取输入视频
videoReader = VideoReader('G:\CA3_rawdata\CA3_2p\B\2024-05-30-1464_merged.mp4');

% 创建输出视频写入对象
videoWriter = VideoWriter('G:\CA3_rawdata\CA3_2p\B\output_video.mp4', 'MPEG-4');
open(videoWriter);

% 设定要保存的帧数
numFrames = min(200, floor(videoReader.FrameRate * videoReader.Duration)); 

% 逐帧读取并写入新视频
frameCount = 0;
while hasFrame(videoReader) && frameCount < numFrames
    frame = readFrame(videoReader);
    writeVideo(videoWriter, frame);
    frameCount = frameCount + 1;
end

% 关闭视频写入对象
close(videoWriter);
disp('前 200 帧视频已保存为 output_video.mp4');
