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
%   'autosave'       - logical，退出时是否保存完整变量（默认 false）
%   'savepath'       - char，保存路径（默认 'autosave_temp.mat'）
%   'intervals'      - m×2 矩阵，限制处理的时间段或帧段
%   'intervals_mode' - 'auto' (default) | 'time' | 'frame'。'auto' 会尝试自动判别 intervals 是时间还是帧号
%
% 返回:
%   position_re_X, position_re_Y - 修正后的坐标

% ---------------- 参数解析 ----------------
p = inputParser;
addParameter(p, 'autosave', false, @(v)islogical(v) || isnumeric(v));
addParameter(p, 'savepath', 'autosave_temp.mat', @ischar);
addParameter(p, 'intervals', [], @(v) isempty(v) || (isnumeric(v) && size(v,2)==2));
addParameter(p, 'intervals_mode', 'auto', @(s)ischar(s) && ismember(lower(s), {'auto','time','frame'}));
parse(p, varargin{:});
autosave = logical(p.Results.autosave);
savepath = p.Results.savepath;
intervals = p.Results.intervals;
intervals_mode = lower(p.Results.intervals_mode);

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

% ---------- 处理 intervals（支持 time 或 frame，两种模式，默认 auto） ----------
% 为了让回调函数访问，modeUsed 与 intervalStarts 提升至此作用域
modeUsed = '';
intervalStarts = []; % 将保存每个 interval 的首帧索引（整数帧号）
origNanIdx = origNanAll; % 默认：所有 NaN
numNaNFiltered = numel(origNanIdx);

if ~isempty(intervals)
    intervals = reshape(intervals, [], 2);
    modeUsed = intervals_mode;
    if strcmp(modeUsed, 'auto')
        allInts = all(mod(intervals(:),1) == 0);
        allWithinFrames = all(intervals(:) >= 1 & intervals(:) <= nFrames);
        if allInts && allWithinFrames
            modeUsed = 'frame';
        else
            modeUsed = 'time';
        end
    end

    switch modeUsed
        case 'time'
            inInterval = false(nFrames,1);
            for k = 1:size(intervals,1)
                t1 = min(intervals(k,:));
                t2 = max(intervals(k,:));
                inInterval = inInterval | (times >= t1 & times <= t2);
            end
            origNanIdx = origNanAll(inInterval(origNanAll));
            fprintf('Intervals interpreted as TIME ranges (s).\n');

            % 计算每个 interval 的首帧（第一个 times >= t1）
            for k = 1:size(intervals,1)
                t1 = min(intervals(k,:));
                % 找到第一个满足 times >= t1 的帧索引
                idx = find(times >= t1, 1, 'first');
                if ~isempty(idx)
                    intervalStarts(end+1,1) = idx; %#ok<AGROW>
                end
            end

        case 'frame'
            % 视作帧索引区间（round 并剪裁）
            inInterval = false(nFrames,1);
            for k = 1:size(intervals,1)
                f1 = max(1, min(nFrames, round(min(intervals(k,:)))));
                f2 = max(1, min(nFrames, round(max(intervals(k,:)))));
                inInterval(f1:f2) = true;
                intervalStarts(end+1,1) = f1; %#ok<AGROW>
            end
            origNanIdx = origNanAll(inInterval(origNanAll));
            fprintf('Intervals interpreted as FRAME index ranges.\n');

        otherwise
            error('Unknown intervals_mode: %s', intervals_mode);
    end

    % 去重并排序 intervalStarts，确保在 1..nFrames
    intervalStarts = unique(intervalStarts(:)');
    intervalStarts = intervalStarts(intervalStarts >= 1 & intervalStarts <= nFrames);
    numNaNFiltered = numel(origNanIdx);
end

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
                % 若没有任何 NaN（origNanAll 为空）且有 intervalStarts，则在 intervalStarts 之间跳转（上一个）
                if isempty(origNanAll) && ~isempty(intervalStarts)
                    prevIdxs = intervalStarts(intervalStarts < currentFrame);
                    if ~isempty(prevIdxs)
                        displayFrame(prevIdxs(end));
                    else
                        % 若无更早的 interval 开始帧，则跳到最后一个 interval 开始（循环）
                        displayFrame(intervalStarts(end));
                    end
                else
                    if currentFrame>1, displayFrame(currentFrame-1); end
                end
            case 'downarrow'
                % 若没有任何 NaN（origNanAll 为空）且有 intervalStarts，则在 intervalStarts 之间跳转（下一个）
                if isempty(origNanAll) && ~isempty(intervalStarts)
                    nextIdxs = intervalStarts(intervalStarts > currentFrame);
                    if ~isempty(nextIdxs)
                        displayFrame(nextIdxs(1));
                    else
                        % 若无更晚的 interval 开始帧，则跳到第一个 interval 开始（循环）
                        displayFrame(intervalStarts(1));
                    end
                else
                    if currentFrame<nFrames, displayFrame(currentFrame+1); end
                end
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
