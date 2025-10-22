function [X_clean, Y_clean, position_timelite, BW1, info] = process_mouse_tracking_oneROI(Path, animal, rec_day, varargin)
% process_mouse_tracking_oneROI 读取视频和轨迹，手动绘制单个活动区域并清理错误点
%
% USAGE:
%   [X_clean, Y_clean, position_timelite, BW1, info] = ...
%       process_mouse_tracking_oneROI(Path, animal, rec_day)
%
% NAME-VALUE OPTIONS:
%   'ConfThresh'    - 置信度阈值（默认 0.8）
%   'DistThresh'    - 相邻点距离阈值（像素，默认 40）
%   'StartFrame'    - 从第几个样本（轨迹行）开始处理（默认 1）
%
% OUTPUTS:
%   X_clean, Y_clean   - 清理后的轨迹（像素，已恢复到原始视频坐标）
%   position_timelite  - 时间向量（秒）对应每个轨迹样本（原始全部）
%   BW1                - cell 数组，包含单个二值掩码 BW1{1}
%   info               - 结构体，包含中间信息（framerate, raw X/Y, poss_nose 等）

%% parse inputs
p = inputParser;
addRequired(p,'Path',@ischar);
addRequired(p,'animal',@ischar);
addRequired(p,'rec_day',@ischar);
addParameter(p,'ConfThresh',0.8,@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,'DistThresh',40,@(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'StartFrame',1,@(x) isnumeric(x) && isscalar(x) && x>=1);
parse(p,Path,animal,rec_day,varargin{:});

confThresh = p.Results.ConfThresh;
distThresh = p.Results.DistThresh;
recordedFrameCount = p.Results.StartFrame;

%% locate video file (*.AVI) and CSV (DLC)
videoDir = fullfile(Path, animal, rec_day, 'video_track');
if ~isfolder(videoDir)
    error('Video track folder not found: %s', videoDir);
end

mp4_file = dir(fullfile(videoDir, '*.AVI'));
if isempty(mp4_file)
    error('No .AVI file found in %s', videoDir);
end
mp4_file = mp4_file(1);

v = VideoReader(fullfile(mp4_file.folder, mp4_file.name));
% MATLAB 的属性名是 FrameRate（注意大小写）
framerate = v.FrameRate;

% 查找 CSV（假定只有一个）
path_file = dir(fullfile(videoDir, '*.csv'));
if isempty(path_file)
    error('No tracking CSV found in %s', videoDir);
end
path_file = path_file(1);
data_path = readtable(fullfile(path_file.folder, path_file.name));

% 时间向量（按行）
position_timelite = (data_path.scorer + 1) / framerate;

% 假定 CSV 列格式与你原脚本一致
poss_nose = table2array(data_path(:,4));            % 第4列置信度
X_raw = table2array(data_path(:,2)) + 100;         % 加100对应 padding
Y_raw = table2array(data_path(:,3)) + 100;

% 读取并填充第一帧（上下左右各100像素）
firstFrame = readFrame(v);
[rows, cols, channels] = size(firstFrame);
newRows = rows + 200;
newCols = cols + 200;
paddedFrame = uint8(255 * ones(newRows, newCols, channels)); % 白色背景
paddedFrame(101:100+rows, 101:100+cols, :) = firstFrame;

%% load or create single BW mask
matfile_name = fullfile(videoDir, 'grab_picture.mat');
jpgfile_name = fullfile(videoDir, 'grab_picture.jpg');

if ~exist(matfile_name,'file')
    % 显示填充后的第一帧并让用户绘制单个ROI
    hFig = figure('Name','Draw single ROI - press Enter when done','NumberTitle','off');
    imshow(paddedFrame);
    hold on;
    scatter(X_raw, Y_raw, 6, 'filled'); % 显示轨迹点参考
    title('Draw a single ROI with mouse (roipoly). Press Enter to finish the polygon.');
    BW = roipoly; % 用户通过鼠标绘制多边形并回车确认
    if isempty(BW)
        warning('No ROI drawn. Exiting and returning empty ROI.');
        BW1 = {};
    else
        BW1 = {BW};
        boundary = bwboundaries(BW);
        if ~isempty(boundary)
            plot(boundary{1}(:,2), boundary{1}(:,1), 'LineWidth', 2);
            mark_x = boundary{1}(1,2);
            mark_y = boundary{1}(1,1);
            text(mark_x, mark_y, '0', 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    % 保存图片和 mat
    try
        saveas(hFig, jpgfile_name, 'jpeg');
        save(matfile_name, 'BW1', 'recordedFrameCount', '-mat');
    catch ME
        warning('Failed to save grab_picture files: %s', ME.message);
    end
    hold off;
    if ishandle(hFig), close(hFig); end
else
    s = load(matfile_name, 'BW1', 'recordedFrameCount');
    if isfield(s, 'BW1')
        BW1 = s.BW1;
    else
        error('grab_picture.mat exists but does not contain BW1');
    end
    if isfield(s, 'recordedFrameCount')
        recordedFrameCount = s.recordedFrameCount;
    end
end

%% 截取轨迹起点并初始化
X_filter = X_raw(recordedFrameCount:end);
Y_filter = Y_raw(recordedFrameCount:end);
poss_filter = poss_nose(recordedFrameCount:end);
time_filter = position_timelite(recordedFrameCount:end);

X_clean = X_filter;
Y_clean = Y_filter;

%% 1) 区域外点设为 NaN （使用第一个 ROI，即 BW1{1}）
if ~isempty(BW1)
    mask = BW1{1};
    mask_rows = size(mask,1);
    mask_cols = size(mask,2);
    idx_valid = ~isnan(X_clean) & ~isnan(Y_clean);
    Xv = X_clean(idx_valid);
    Yv = Y_clean(idx_valid);
    xi = fix(Yv + 1); % Y -> 行索引
    yi = fix(Xv + 1); % X -> 列索引
    out_of_bounds = xi < 1 | xi > mask_rows | yi < 1 | yi > mask_cols;
    valid_linear = true(size(xi));
    valid_linear(out_of_bounds) = false;
    lin_idx = nan(size(xi));
    lin_idx(valid_linear) = sub2ind(size(mask), xi(valid_linear), yi(valid_linear));
    mask_values = false(size(xi));
    mask_values(valid_linear) = mask(lin_idx(valid_linear));
    invalid_indices_global = find(idx_valid);
    invalid_mask_positions = invalid_indices_global(~mask_values);
    X_clean(invalid_mask_positions) = NaN;
    Y_clean(invalid_mask_positions) = NaN;
end

%% 2) 置信度过滤
lowconf = poss_filter < confThresh;
X_clean(lowconf) = NaN;
Y_clean(lowconf) = NaN;

%% 3) 相邻点距离过滤：若前后距离 > distThresh，则后者设为 NaN
n = numel(X_clean);
for i = 2:n
    if isnan(X_clean(i-1)) || isnan(Y_clean(i-1)) || isnan(X_clean(i)) || isnan(Y_clean(i))
        continue;
    end
    d = hypot( X_clean(i) - X_clean(i-1), Y_clean(i) - Y_clean(i-1) );
    if d > distThresh
        X_clean(i) = NaN;
        Y_clean(i) = NaN;
    end
end

%% 恢复到原始视频坐标（去掉 +100 padding）
position_X = X_clean - 100;
position_Y = Y_clean - 100;

X_clean = position_X;
Y_clean = position_Y;

%% info 输出
info = struct();
info.framerate = framerate;
info.rawX = X_raw - 100;
info.rawY = Y_raw - 100;
info.poss_nose = poss_nose;
info.time_all = position_timelite;
info.time_filtered = time_filter;

end
