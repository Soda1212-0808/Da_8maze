function display_next_frame_on_scroll(videoFile)
    % 检查视频文件是否存在
    if ~exist(videoFile, 'file')
        error('The specified video file does not exist.');
    end
    
    % 打开视频文件
    videoReader = VideoReader(videoFile);
    
    % 获取视频的帧率
    frameRate = videoReader.FrameRate;
    
    % 初始化帧计数器
    frameCount = 0;
    
    % 创建图形窗口
    hFig = figure;
    set(hFig, 'UserData', struct('videoReader', videoReader, 'frameCount', frameCount, 'frameRate', frameRate));
    
    % 设置滚轮事件和键盘按键回调函数
    set(hFig, 'WindowScrollWheelFcn', @(src, event) onScroll(src, event));
    set(hFig, 'KeyPressFcn', @(src, event) onKeyPress(src, event));
    
    % 显示第一帧
    if hasFrame(videoReader)
        frame = readFrame(videoReader);
        frameCount = frameCount + 1;
        set(hFig, 'UserData', struct('videoReader', videoReader, 'frameCount', frameCount, 'frameRate', frameRate));
        imshow(frame, 'Parent', gca);
        % 显示帧数
        displayFrameCount(frameCount);
    end
end

function onScroll(src, event)
    data = get(src, 'UserData');
    videoReader = data.videoReader;
    frameCount = data.frameCount;
    frameRate = data.frameRate;
    
    % 检查滚轮方向
    if event.VerticalScrollCount < 0
        % 滚轮向上滚动，显示下一帧
        frameCount = frameCount + 1;
    elseif event.VerticalScrollCount > 0
        % 滚轮向下滚动，显示前一帧
        frameCount = frameCount - 1;
    end
    
    % 更新帧计数器
    if frameCount < 1
        frameCount = 1;
    elseif frameCount > floor(videoReader.Duration * frameRate)
        frameCount = floor(videoReader.Duration * frameRate);
    end
    
    % 跳转到新时间点并读取帧
    videoReader.CurrentTime = (frameCount - 1) / frameRate;
    frame = readFrame(videoReader);
    
    % 更新帧计数器
    data.frameCount = frameCount;
    set(src, 'UserData', data);
    
    % 显示当前帧
    imshow(frame, 'Parent', gca);
    % 显示帧数
    displayFrameCount(frameCount);
end

function onKeyPress(src, event)
    if strcmp(event.Key, 'return')
        % 获取当前帧计数器
        data = get(src, 'UserData');
        frameCount = data.frameCount;
        
        % 记录当前帧数
        disp(['Recorded frame number: ' num2str(frameCount)]);
        
        % 保存当前帧数到工作空间变量
        assignin('base', 'recordedFrameCount', frameCount);
        
        % 关闭图像窗口
        close(src);
    end
end

function displayFrameCount(frameCount)
    % 在图像右上角显示帧数
    text(10, 10, ['Frame: ' num2str(frameCount)], 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'black', 'Margin', 2);
end
