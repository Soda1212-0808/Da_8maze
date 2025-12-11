% 推荐：不要用 clear all，改成有选择地清理
clear all
close all

animal = 'DCA3-9';
% rec_day = '2021-06-13';
recordings=ds.find_recordings(animal);
rec_day = recordings(2).day;

Path = ds.locations.server_data_path;

% 参数（把可配置项写到顶部，方便调整）
confThresh = 0.8;
distThresh = 40;
startFrame = 1;
videoTrackDir = fullfile(Path, animal, rec_day, 'video_track');

% 检查目录
if ~isfolder(videoTrackDir)
    error('视频轨迹目录不存在：%s', videoTrackDir);
end

% 读取并处理轨迹（process_mouse_tracking_oneROI 或 process_mouse_tracking）
try
    [position_X, position_Y, position_timelite, BW, info] = ...
        process_mouse_tracking(Path, animal, rec_day, ...
            'ConfThresh', confThresh, 'DistThresh', distThresh, 'StartFrame', startFrame);
catch ME
    error('调用 process_mouse_tracking 出错：%s\n%s', ME.message, ME.getReport());
end

% 保存清理后的轨迹（用 -v7.3 当数据可能较大）
outfile = fullfile(videoTrackDir, 'mice_position.mat');
try
    save(outfile, 'position_X', 'position_Y', 'position_timelite', '-v7.3');
catch ME
    warning('保存 mice_position.mat 失败：%s', ME.message);
    save(outfile, 'position_X', 'position_Y', 'position_timelite'); % 尝试普通方式
end

% 手动修正保存（改进后）
% 找到 behavior csv
behav_csv = dir(fullfile(Path, animal, rec_day, 'behavior', '*.csv'));
if isempty(behav_csv)
    warning('未找到行为 CSV 文件，跳过 intervals 加载');
    intervals = [];
else
    try
        data_event = readtable(fullfile(behav_csv(1).folder, behav_csv(1).name));
        data_evt_arr = table2array(data_event);
        % 只在有足够列时构造 intervals
        if size(data_evt_arr,2) >= 9
            intervals = [data_evt_arr(:,3) data_evt_arr(:,5) ; data_evt_arr(:,7) data_evt_arr(:,9)];
        else
            intervals = [];
            warning('behavior CSV 列数不足，intervals 留空');
        end
    catch ME
        warning('读取 behavior CSV 失败：%s', ME.message);
        intervals = [];
    end
end

% 准备 mice_position_re 文件，如果存在就加载，否则创建
pos_re_file = fullfile(videoTrackDir, 'mice_position_re.mat');
if isfile(pos_re_file)
    S = load(pos_re_file);
    if isfield(S, 'position_re_X') && isfield(S, 'position_re_Y')
        position_re_X = S.position_re_X;
        position_re_Y = S.position_re_Y;
    else
        position_re_X = position_X;
        position_re_Y = position_Y;
    end
else
    position_re_X = position_X;
    position_re_Y = position_Y;
    save(pos_re_file, 'position_re_X', 'position_re_Y', 'position_timelite', '-v7.3');
end

%%

% 寻找视频文件（缓存）
video_file = dir(fullfile(videoTrackDir, '*.AVI'));
if isempty(video_file)
    video_file = dir(fullfile(videoTrackDir, '*.mp4')); % 备选
end
if isempty(video_file)
    error('未找到视频文件（.AVI or .mp4）在 %s', videoTrackDir);
end
video_fullpath = fullfile(video_file(1).folder, video_file(1).name);

% 调用交互修正函数（manual_fix_nan_frames）
[position_re_X, position_re_Y] = manual_fix_nan_frames(video_fullpath, ...
    position_re_X, position_re_Y, position_timelite, ...
    'intervals', intervals, 'autosave', true, 'savepath', pos_re_file);
