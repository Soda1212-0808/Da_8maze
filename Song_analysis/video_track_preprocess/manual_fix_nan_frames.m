function [position_re_X, position_re_Y] = manual_fix_nan_frames(videoFile, x, y, position_timelite, varargin)
% manual_fix_nan_frames8
% 手动修正视频中检测失败（NaN）的帧坐标（必含时间轴，自动保存带 timelist）
%
% 必要输入:
%   videoFile          - 视频文件路径（string）
%   x, y               - n×1 向量（目标点坐标）
%   position_timeline  - n×1 向量，每帧对应的时间（单位：秒）
%
% 可选 Name-Value 参数:
%   'autosave'  - logical，退出时是否保存完整变量（默认 false）
%   'savepath'  - char，保存路径（默认 'autosave_temp.mat'）
%   'intervals' - m×2 矩阵，限制处理的时间段 [t_start, t_end]
%
% 返回:
%   x_new, y_new - 修正后的坐标
%
% 示例:
%   % 普通模式（处理所有 NaN 帧）
%   [x2,y2] = manual_fix_nan_frames8('video.mp4', x, y, tlist);
%
%   % 限定时间段 + 自动保存
%   [x2,y2] = manual_fix_nan_frames8('video.mp4', x, y, tlist, ...
%       'intervals', [5 10; 20 25], 'autosave', true, 'savepath', 'progress.mat');

% ---------------- 参数解析 ----------------
p = inputParser;
addParameter(p, 'autosave', false, @(v)islogical(v) || isnumeric(v));
addParameter(p, 'savepath', 'autosave_temp.mat', @ischar);
addParameter(p, 'intervals', [], @(v) isempty(v) || (isnumeric(v) && size(v,2)==2));
parse(p, varargin{:});
autosave = logical(p.Results.autosave);
savepath = p.Results.savepath;
intervals = p.Results.intervals;

% ---------- 输入检查 ----------
if numel(x) ~= numel(y)
    error('x and y must have same length.');
end
nFrames = numel(x);
if numel(position_timelite) ~= nFrames
    error('position_timelist length (%d) must match x,y length (%d).', ...
        numel(position_timelite), nFrames);
end
position_re_X = x(:);
position_re_Y = y(:);
times = position_timelite(:);

% ---------- 读取视频 ----------
v = VideoReader(videoFile);
frameRate = v.FrameRate;
try
    totalFramesVideo = v.NumFrames; %#ok<NASGU>
catch
    totalFramesVideo = floor(v.Duration * frameRate);
end

% ---------- 计算 NaN 帧集合 ----------
origNanAll = find(isnan(x) | isnan(y));
if isempty(origNanAll)
    warning('No NaN frames found in x or y.');
end

% 若指定 intervals，则筛选时间段内的 NaN 帧
if ~isempty(intervals)
    inInterval = false(nFrames,1);
    for k = 1:size(intervals,1)
        t1 = min(intervals(k,:));
        t2 = max(intervals(k,:));
        inInterval = inInterval | (times >= t1 & times <= t2);
    end
    origNanIdx = origNanAll(inInterval(origNanAll));
else
    origNanIdx = origNanAll;
end
numNaNFiltered = numel(origNanIdx);

fprintf('Total frames: %d | NaN frames: %d | Filtered NaN (in intervals): %d\n', ...
    nFrames, numel(origNanAll), numNaNFiltered);

% ---------- GUI 初始化 ----------
currentFrame = 1;
hFig = figure('Name','Manual fix NaN frames (FullScreen)','NumberTitle','off',...
    'KeyPressFcn',@keyPressCallback);
set(hFig,'Units','normalized','OuterPosition',[0 0 1 1]);
hAx = axes('Parent',hFig,'Position',[0.05 0.08 0.9 0.85]);
axis(hAx,'image');
hold(hAx,'on');

displayFrame(currentFrame);
uiwait(hFig);

% ---------- 退出时保存 ----------
if autosave
    try
        save(savepath, 'position_re_X', 'position_re_Y', 'position_timelite');
        fprintf('Final save: %s\n', savepath);
    catch ME
        warning('Final save failed: %s', ME.message);
    end
end
return

%% ============ 内部函数 ============

    function displayFrame(frameNum)
        frameNum = max(1,min(nFrames,frameNum));
        currentFrame = frameNum;
        cla(hAx);

        try
            frame = read(v, frameNum);
        catch
            v.CurrentTime = max(0,(frameNum-1)/frameRate);
            frame = readFrame(v);
        end
        imh = imshow(frame,'Parent',hAx);
        set(imh,'ButtonDownFcn',@imageClickCallback);

        % 绘制已有标记
        if ~isnan(position_re_X(frameNum)) && ~isnan(position_re_Y(frameNum))
            plot(hAx, position_re_X(frameNum), position_re_Y(frameNum), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
        end

        % 构建标题
        isNanFiltered = ismember(frameNum, origNanIdx);
        if isNanFiltered
            idx = find(origNanIdx == frameNum);
            nanInfo = sprintf('NaN frame #%d / %d (filtered)', idx, numNaNFiltered);
        else
            nanInfo = sprintf('Non-NaN frame (filtered total %d)', numNaNFiltered);
        end
        tcur = times(frameNum);
        titleStr = sprintf(['Frame %d / %d | %s\n',...
            'Time = %.3f s | x = %.2f, y = %.2f\n',...
            '(↑/↓: prev/next, ←/→: filtered-NaN jump, click=set, q=quit)'], ...
            frameNum, nFrames, nanInfo, tcur, position_re_X(frameNum), position_re_Y(frameNum));
        title(hAx, titleStr,'Interpreter','none','FontSize',11,'FontWeight','bold');
        drawnow;
    end

    function imageClickCallback(~,~)
        cp = get(hAx,'CurrentPoint');
        xg = cp(1,1); yg = cp(1,2);
        position_re_X(currentFrame) = xg;
        position_re_Y(currentFrame) = yg;
        fprintf('Frame %d marked (%.2f, %.2f). Auto-saving...\n', currentFrame, xg, yg);
        try
            save(savepath, 'position_re_X', 'position_re_Y', 'position_timelite');
        catch ME
            warning('Auto-save failed: %s', ME.message);
        end
        displayFrame(currentFrame);
    end

    function keyPressCallback(~,event)
        switch event.Key
            case 'uparrow'
                if currentFrame>1, displayFrame(currentFrame-1); end
            case 'downarrow'
                if currentFrame<nFrames, displayFrame(currentFrame+1); end
            case 'leftarrow'
                leftIdx = origNanIdx(origNanIdx<currentFrame);
                if ~isempty(leftIdx), displayFrame(max(leftIdx)); end
            case 'rightarrow'
                rightIdx = origNanIdx(origNanIdx>currentFrame);
                if ~isempty(rightIdx), displayFrame(min(rightIdx)); end
            case 'q'
                fprintf('Quit requested. Closing window...\n');
                uiresume(hFig);
                close(hFig);
        end
    end
end
