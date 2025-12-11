function [position_re_X, position_re_Y] = manual_fix_nan_frames2(videoFile, x, y, position_timelite, varargin)
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
%
% 说明:
%   - intervals 有两种可能的语义:
%       时间段（秒）： [t_start, t_end]
%       帧段（帧索引）： [f_start, f_end]，帧索引从 1 开始到 nFrames
%   - 若 intervals_mode = 'auto'，将按下列规则判别：
%       * 若 intervals 中所有元素均为整数且都位于 [1, nFrames]，则认为是帧段（frame）
%       * 否则认为是时间段（time）
%
% 示例:
%   [x2,y2] = manual_fix_nan_frames8('video.mp4', x, y, tlist);
%   [x2,y2] = manual_fix_nan_frames8('video.mp4', x, y, tlist, ...
%       'intervals', [5 10; 20 25], 'autosave', true, 'savepath', 'progress.mat');
%   % 用帧号区间（强制 frame）
%   [x2,y2] = manual_fix_nan_frames8('video.mp4', x, y, tlist, ...
%       'intervals', [100 200; 500 550], 'intervals_mode', 'frame');

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
origNanIdx = origNanAll; % 默认全部 NaN
numNaNFiltered = numel(origNanIdx);

%% ADDED: 初始化 intervalStarts（用于无 NaN 时在 intervals 之间跳转）
intervalStarts = [];  % 将保存每个 interval 的起始帧（frame index）

if ~isempty(intervals)
    % 确保 intervals 是 m x 2
    intervals = reshape(intervals, [], 2);
    % 判别模式
    modeUsed = intervals_mode;
    if strcmp(modeUsed, 'auto')
        % 自动检测：如果 intervals 所有元素都是整数并且都在 [1, nFrames]，则当作帧索引
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
            % 原来基于时间的筛选
            inInterval = false(nFrames,1);
            for k = 1:size(intervals,1)
                t1 = min(intervals(k,:));
                t2 = max(intervals(k,:));
                inInterval = inInterval | (times >= t1 & times <= t2);
            end
            origNanIdx = origNanAll(inInterval(origNanAll));
            fprintf('Intervals interpreted as TIME ranges (s).\n');
        case 'frame'
            % 将 intervals 视作帧索引区间（可能是浮点数，但会 round）
            inInterval = false(nFrames,1);
            for k = 1:size(intervals,1)
                % 四舍五入并限制在 [1, nFrames]
                f1 = max(1, min(nFrames, round(min(intervals(k,:)))));
                f2 = max(1, min(nFrames, round(max(intervals(k,:)))));
                inInterval(f1:f2) = true;
            end
            origNanIdx = origNanAll(inInterval(origNanAll));
            fprintf('Intervals interpreted as FRAME index ranges.\n');
        otherwise
            error('Unknown intervals_mode: %s', intervals_mode);
    end
    numNaNFiltered = numel(origNanIdx);

    %% ADDED: 根据 modeUsed 填充 intervalStarts（每个 interval 的第一帧）
    try
        if exist('modeUsed','var') && strcmp(modeUsed,'frame')
            % intervals 每行为帧号区间，取每行的起始帧（四舍五入并裁剪）
            for k = 1:size(intervals,1)
                f1 = round(min(intervals(k,:)));
                f1 = max(1, min(nFrames, f1));
                f2 = round(max(intervals(k,:))); %#ok<NASGU>
                if f1 <= nFrames
                    intervalStarts(end+1,1) = f1; %#ok<AGROW>
                end
            end
        else
            % time 模式：取每个时间区间的第一个匹配帧（times >= t1）
            for k = 1:size(intervals,1)
                t1 = min(intervals(k,:));
                t2 = max(intervals(k,:));
                idxStart = find(times >= t1, 1, 'first');
                if isempty(idxStart)
                    idxStart = find(times <= t2, 1, 'first');
                end
                if ~isempty(idxStart)
                    intervalStarts(end+1,1) = idxStart; %#ok<AGROW>
                end
            end
        end
    catch
        % 如果在计算 intervalStarts 时出错，保留为空（不影响其他逻辑）
        intervalStarts = [];
    end

    % 规范化 intervalStarts：去重、升序、裁剪为合法帧索引
    if ~isempty(intervalStarts)
        intervalStarts = floor(intervalStarts(:));
        intervalStarts = max(1, min(nFrames, intervalStarts));
        intervalStarts = unique(intervalStarts, 'sorted')';
    end
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
        % DEBUG 增强版 displayFrame：打印、捕获错误、强制 drawnow/pause
        frameNum = max(1,min(nFrames,frameNum));
        fprintf('DEBUG: entering displayFrame(%d) (currentFrame before=%d)\n', frameNum, currentFrame);
        currentFrame = frameNum;
        cla(hAx);

        try
            % 尝试直接用 read（若失败走 catch 分支）
            try
                frame = read(v, frameNum);
            catch readErr
                % 如果 read 失败，尝试通过设置 CurrentTime 并 readFrame
                fprintf('DEBUG: read(v,%d) failed: %s. Trying seek via CurrentTime/readFrame...\n', frameNum, readErr.message);
                v.CurrentTime = max(0,(frameNum-1)/frameRate);
                frame = readFrame(v);
            end
        catch ME
            % 如果仍然失败，打印错误并返回（不 crash）
            fprintf('ERROR: failed to read frame %d: %s\n', frameNum, ME.message);
            % 显示空白或提示文本
            cla(hAx);
            text(0.5,0.5,sprintf('Failed to read frame %d', frameNum),'Units','normalized','HorizontalAlignment','center','Parent',hAx);
            drawnow;
            return;
        end

        % 显示帧（并检查 frame 是否为空）
        if isempty(frame)
            fprintf('DEBUG: read returned empty frame for %d\n', frameNum);
            cla(hAx);
            text(0.5,0.5,sprintf('Empty frame %d', frameNum),'Units','normalized','HorizontalAlignment','center','Parent',hAx);
            drawnow;
            return;
        end

        imh = imshow(frame,'Parent',hAx);
        set(imh,'ButtonDownFcn',@imageClickCallback);

        % 绘制已有标记
        if ~isnan(position_re_X(frameNum)) && ~isnan(position_re_Y(frameNum))
            hold(hAx,'on');
            plot(hAx, position_re_X(frameNum), position_re_Y(frameNum), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
            hold(hAx,'off');
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

        % 强制刷新并短暂停顿，确保 GUI 更新
        drawnow;
        pause(0.01);
        fprintf('DEBUG: displayed frame %d successfully\n', frameNum);
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
        % 更稳健的 keyPressCallback（仅替换此处）
        fprintf('KEYPRESS: %s | currentFrame=%d | numNaN=%d | intervalStarts=%s\n', ...
            event.Key, currentFrame, numel(origNanAll), mat2str(intervalStarts));

        switch event.Key
            case 'uparrow'
                if currentFrame>1, displayFrame(currentFrame-1); end
            case 'downarrow'
                if currentFrame<nFrames, displayFrame(currentFrame+1); end
            case 'leftarrow'
                % 判断应使用哪个集合来查找上一个位置
                if ~isempty(origNanAll)
                    % 优先使用经过 intervals 筛选的 NaN（origNanIdx）
                    searchSet = origNanIdx;
                    % 若 filtered 为空，则优先使用 intervalStarts（若存在），否则回退到 origNanAll
                    if isempty(searchSet)
                        if ~isempty(intervalStarts)
                            % treat intervalStarts as possible jump targets
                            idx = find(intervalStarts < currentFrame, 1, 'last');
                            if ~isempty(idx)
                                displayFrame(intervalStarts(idx));
                            else
                                % 没有更早的 interval start -> 跳到第一个 interval start
                                displayFrame(intervalStarts(1));
                            end
                            return;
                        else
                            searchSet = origNanAll;
                        end
                    end
                    % 用 searchSet 做 NaN 跳转（若仍为空则什么也不做）
                    if ~isempty(searchSet)
                        leftIdx = searchSet(searchSet < currentFrame);
                        if ~isempty(leftIdx)
                            displayFrame(max(leftIdx));
                        end
                    end
                else
                    % 没有 NaN，按 intervalStarts 或逐帧退回
                    if ~isempty(intervalStarts)
                        idx = find(intervalStarts < currentFrame, 1, 'last');
                        if ~isempty(idx), displayFrame(intervalStarts(idx));
                        else
                            if currentFrame>1, displayFrame(currentFrame-1); end
                        end
                    else
                        if currentFrame>1, displayFrame(currentFrame-1); end
                    end
                end

            case 'rightarrow'
                % 判断应使用哪个集合来查找下一个位置
                if ~isempty(origNanAll)
                    % 优先使用经过 intervals 筛选的 NaN（origNanIdx）
                    searchSet = origNanIdx;
                    % 若 filtered 为空，则优先使用 intervalStarts（若存在），否则回退到 origNanAll
                    if isempty(searchSet)
                        if ~isempty(intervalStarts)
                            idx = find(intervalStarts > currentFrame, 1, 'first');
                            if ~isempty(idx)
                                displayFrame(intervalStarts(idx));
                            else
                                % 没有更晚的 interval start -> 跳到最后一个 interval start
                                displayFrame(intervalStarts(end));
                            end
                            return;
                        else
                            searchSet = origNanAll;
                        end
                    end
                    % 用 searchSet 做 NaN 跳转（若仍为空则什么也不做）
                    if ~isempty(searchSet)
                        rightIdx = searchSet(searchSet > currentFrame);
                        if ~isempty(rightIdx)
                            displayFrame(min(rightIdx));
                        end
                    end
                else
                    % 没有 NaN，按 intervalStarts 或逐帧前进
                    if ~isempty(intervalStarts)
                        idx = find(intervalStarts > currentFrame, 1, 'first');
                        if ~isempty(idx), displayFrame(intervalStarts(idx));
                        else
                            if currentFrame<nFrames, displayFrame(currentFrame+1); end
                        end
                    else
                        if currentFrame<nFrames, displayFrame(currentFrame+1); end
                    end
                end

            case 'q'
                fprintf('Quit requested. Closing window...\n');
                uiresume(hFig);
                close(hFig);
        end
    end


end
