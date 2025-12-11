function results = detect_place_cell_v3(positionX, positionY, position_timelite, spike_times, inIntervals, bin_size, varargin)
% DETECT_PLACE_CELL_V2 (optimized, heavily commented, bilingual)
% Detect place fields and test whether a neuron is a place cell.
% 支持 place field 检测及显著性检验
%
% Modes / 模式:
%   'spike'  : spike_times is a vector of event timestamps (seconds).  （事件时间戳：秒）
%   'ca'     : spike_times is actually a Tx1 activity vector aligned to position_timelite
%              (e.g. deconvolved calcium amplitude per frame). 在帧对齐下的去卷积活动值（按帧累计）
%
% Major optimizations / 主要优化点:
%   - single-call interpolation for spikes instead of per-interval repeated interp1
%     （对 spikes 只做一次插值，而不是对每个区段重复插值）
%   - discretize + accumarray for per-bin accumulation (fast)
%     （使用 discretize + accumarray 来做箱内累计，效率高）
%   - vectorized occupancy accumulation（向量化计算占用时间）
%   - optional parfor for shuffles (useParpool)（对 shuffle 循环可选并行）
%   - optional GPU separable smoothing (useGPU)（可选 GPU 可分离卷积平滑）
%   - optional single precision (useSingle) to reduce memory/use GPU faster
%     （可选 single 精度，节省内存并加速 GPU）
%
% USAGE:
%   results = detect_place_cell_v2(posX,posY,posTime,spikesOrActivity,inIntervals,bin_size, Name,Value,...)
%   （用法示例）
%
% See inline comments for detailed bilingual explanation of each block.
% 详见代码内注释（中英双语）。
% -------------------------------------------------------------------------

% -----------------------------
% Parse name-value inputs
% 解析名-值参数
% -----------------------------
p = inputParser;
addParameter(p,'smooth',false,@(v)islogical(v)||isnumeric(v));                 % whether to smooth / 是否平滑
addParameter(p,'smooth_bin',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);             % smoothing scale in bins / 平滑尺度（箱数）
addParameter(p,'bin_pass',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);              % minimum sample-count per bin / 最小样本数阈值
addParameter(p,'min_t',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=0));% alternative to set bin_pass / 可用 min_t 设置 bin_pass
addParameter(p,'frame_rate',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));% required for min_t / 若用 min_t 则需 frame_rate
addParameter(p,'nShuffles',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);         % number of shuffles / 置换次数
addParameter(p,'info_thresh',0.5,@(x)isnumeric(x)&&isscalar(x));              % info threshold bits/spike / 信息量阈值
addParameter(p,'p_thresh',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);      % p-value threshold / p 值阈值
addParameter(p,'min_field_bins',3,@(x)isnumeric(x)&&isscalar(x)&&x>0);        % minimal field size in bins / 最小场大小（箱数）
addParameter(p,'peak_rate_min',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);        % minimal peak rate for field / 峰率阈值
addParameter(p,'field_threshold_fraction',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1); % field threshold fraction / 场阈值比例
addParameter(p,'verbose',true,@(v)islogical(v)||isnumeric(v));                % verbosity / 是否打印
addParameter(p,'data_type','spike',@(s)ischar(s) || isstring(s));             % 'spike' or 'ca' / 数据类型
addParameter(p,'useParpool',false,@(v)islogical(v)||isnumeric(v));            % use parfor for shuffles / 是否并行
addParameter(p,'useGPU',false,@(v)islogical(v)||isnumeric(v));                % use GPU for smoothing / 是否使用 GPU 平滑
addParameter(p,'useSingle',false,@(v)islogical(v)||isnumeric(v));             % cast to single precision / 是否用 single
parse(p,varargin{:});

% assign options / 赋值本地变量
doSmooth = logical(p.Results.smooth);
smooth_bin = p.Results.smooth_bin;
bin_pass = p.Results.bin_pass;
min_t = p.Results.min_t;
frame_rate = p.Results.frame_rate;
nShuffles = p.Results.nShuffles;
info_thresh = p.Results.info_thresh;
p_thresh = p.Results.p_thresh;
min_field_bins = p.Results.min_field_bins;
peak_rate_min = p.Results.peak_rate_min;
thr_frac = p.Results.field_threshold_fraction;
verbose = logical(p.Results.verbose);
data_type = lower(string(p.Results.data_type));
useParpool = logical(p.Results.useParpool);
useGPU = logical(p.Results.useGPU);
useSingle = logical(p.Results.useSingle);

% compute bin_pass from min_t if provided / 如果提供 min_t 则计算 bin_pass
if ~isempty(min_t)
    if isempty(frame_rate)
        error('frame_rate must be provided when min_t is set. / 当设置 min_t 时必须提供 frame_rate。');
    end
    bin_pass = ceil(min_t * frame_rate);
end

% -----------------------------
% Basic input validation and reshape
% 基本输入检查与 reshape
% -----------------------------
positionX = positionX(:);         % vectorize / 向量化
positionY = positionY(:);
position_timelite = position_timelite(:);
spike_times = spike_times(:);     % may be timestamps or activity vector / 可能是事件时间或 CA 活动向量
T = numel(positionX);
if numel(positionY) ~= T || numel(position_timelite) ~= T
    error('positionX/Y/timelite must have same length T. / 位置与时间向量长度必须相同。');
end
if nargin < 5 || isempty(inIntervals)
    inIntervals = true(T,1);      % default: all samples included / 默认全部样本包含
end
if ~islogical(inIntervals)
    inIntervals = inIntervals ~= 0;
end
if size(inIntervals,1) ~= T
    error('inIntervals must have T rows. / inIntervals 必须有 T 行。');
end

% CA mode validation / CA 模式验证
if data_type == "ca" || data_type == "ca_2p"
    if numel(spike_times) ~= T
        error('For data_type="ca", the 4th input must be a Tx1 activity vector aligned with position_timelite. / CA 模式下，第4个输入必须是与位置时间对齐的 Tx1 活动向量。');
    end
    activity_vec = spike_times(:); % rename / 重命名为 activity_vec
else
    activity_vec = [];
end

% optional cast to single to save memory / 可选使用 single 精度以减少内存
if useSingle
    positionX = single(positionX);
    positionY = single(positionY);
    position_timelite = single(position_timelite);
    if ~isempty(activity_vec), activity_vec = single(activity_vec); end
end

% -----------------------------
% Build spatial bins & occupancy
% 构建空间箱与占用时间
% -----------------------------
valid_pos = ~isnan(positionX) & ~isnan(positionY);   % which samples have valid positions / 有效位置样本
mask_any = any(inIntervals,2);                       % samples belonging to any interval / 属于任何区段的样本
mask_edges = valid_pos & mask_any;
if ~any(mask_edges)
    mask_edges = valid_pos; % fallback / 回退
end

% determine edges with padding so bin centers align well with data / 确定边界并加 pad
xmin = min(positionX(mask_edges)); xmax = max(positionX(mask_edges));
ymin = min(positionY(mask_edges)); ymax = max(positionY(mask_edges));
pad = bin_size;
x_edges = xmin-pad : bin_size : xmax+pad;
y_edges = ymin-pad : bin_size : ymax+pad;
nx = numel(x_edges)-1;
ny = numel(y_edges)-1;

% dt: time per sample (last uses median diff) / 每个位置样本对应时间（最后一项用中位差）
dt = [diff(position_timelite); median(diff(position_timelite))];
if ~any(dt>0)
    dt(:) = 1; % fallback / 回退
end

% include_mask: only those samples with valid pos AND in any interval / 仅包含有效且在区段内的样本
include_mask = valid_pos & mask_any;

% occupancy_time: seconds spent in each bin, computed via discretize + accumarray for speed
occupancy_time = zeros(nx, ny);
if any(include_mask)
    [~,~,xbin_idx] = histcounts(positionX(include_mask), x_edges); % x bin idx for included samples
    [~,~,ybin_idx] = histcounts(positionY(include_mask), y_edges); % y bin idx for included samples
    inds = find(include_mask); % original indices of included samples
    valid_samples = (xbin_idx >= 1) & (ybin_idx >= 1); % discard out-of-edge points
    if any(valid_samples)
        lin_idx = sub2ind([nx, ny], xbin_idx(valid_samples), ybin_idx(valid_samples));
        occ_lin = accumarray(lin_idx, dt(inds(valid_samples)), [nx*ny, 1]); % sum dt per linear index
        occupancy_time = reshape(occ_lin, [nx, ny]);
    else
        occupancy_time = zeros(nx, ny);
    end
else
    occupancy_time = zeros(nx, ny);
end

% sample_counts: number of position samples per bin (used for bin_pass threshold)
if any(include_mask)
    sample_counts = histcounts2(positionX(include_mask), positionY(include_mask), x_edges, y_edges);
else
    sample_counts = zeros(nx, ny);
end

% -----------------------------
% Precompute time intervals from inIntervals (continuous segments)
% 将 inIntervals 转为连续时间区段列表（开始/结束）
% -----------------------------
nCols = size(inIntervals,2);
intervalStarts = [];
intervalEnds = [];
for c = 1:nCols
    col = inIntervals(:,c);
    if ~any(col), continue; end
    d = diff([0; col(:); 0]);
    s = find(d == 1);
    e = find(d == -1) - 1;
    intervalStarts = [intervalStarts; s];
    intervalEnds = [intervalEnds; e];
end
% convert sample indices to times for membership tests
intervalT0 = position_timelite(intervalStarts);
intervalT1 = position_timelite(intervalEnds);

% membership helper (vectorized over intervals) / 判断某些时间是否在任一区间内
is_in_intervals = @(t) any((t(:) >= intervalT0(:)') & (t(:) <= intervalT1(:)'), 2);

% -----------------------------
% Map spikes or CA activity to spatial bins (vectorized)
% 将 spikes 或 CA 活动映射到空间箱（向量化）
% -----------------------------
if data_type == "spike"
    % ---------- Spike mode ----------
    S = numel(spike_times);
    if S == 0
        spk_map = zeros(nx, ny);
    else
        % interpolate spike times once onto position timeline (fast) / 只插值一次
        spk_x_all = interp1(position_timelite, positionX, spike_times, 'linear', NaN);
        spk_y_all = interp1(position_timelite, positionY, spike_times, 'linear', NaN);

        % mask spikes outside inIntervals / 去除不属于任何区段的 spikes
        if ~isempty(intervalT0)
            inMask = false(S,1);
            for k = 1:numel(intervalT0)
                inMask = inMask | (spike_times >= intervalT0(k) & spike_times <= intervalT1(k));
            end
        else
            inMask = false(S,1);
        end
        spk_x_all(~inMask) = NaN;
        spk_y_all(~inMask) = NaN;

        % discretize coordinates then accumulate counts via accumarray (fast)
        xb = discretize(spk_x_all, x_edges);
        yb = discretize(spk_y_all, y_edges);
        valid_sp = ~isnan(xb) & ~isnan(yb);
        if any(valid_sp)
            linIdx = sub2ind([nx, ny], xb(valid_sp), yb(valid_sp));
            counts_lin = accumarray(linIdx, 1, [nx*ny, 1]);
            spk_map = reshape(counts_lin, [nx, ny]);
        else
            spk_map = zeros(nx, ny);
        end
    end
    activity_map_raw = spk_map; % counts (will be converted to rate by dividing occupancy_time)
else
    % ---------- CA mode ----------
    if ~any(include_mask)
        activity_map_raw = zeros(nx, ny);
    else
        inds = find(include_mask);
        xbin_full = nan(T,1);
        ybin_full = nan(T,1);
        xbin_full(inds) = xbin_idx;
        ybin_full(inds) = ybin_idx;
        valid_samples = ~isnan(xbin_full) & ~isnan(ybin_full) & ~isnan(activity_vec);
        if any(valid_samples)
            lin_idx = sub2ind([nx, ny], xbin_full(valid_samples), ybin_full(valid_samples));
            sums = accumarray(lin_idx, double(activity_vec(valid_samples)), [nx*ny, 1]); % sum activity
            activity_map_raw = reshape(sums, [nx, ny]);
        else
            activity_map_raw = zeros(nx, ny);
        end
    end
    spk_map = nan(nx, ny); % undefined in CA mode / CA 模式下 spike_map 设为 NaN
end

% -----------------------------
% Convert raw counts/sums to rate (activity per second)
% 将计数/累加值转为每秒率（activity/sec）
% -----------------------------
rate_map_raw = activity_map_raw ./ occupancy_time; % divide by seconds in each bin
rate_map_raw(~isfinite(rate_map_raw)) = NaN;       % convert Inf/ -Inf -> NaN

% apply bin_pass threshold: mark low-visit bins as zero (optional)
if bin_pass > 0
    low_mask = sample_counts < bin_pass;
    rate_map_raw(low_mask) = 0;
end

% -----------------------------
% Smoothing (separable Gaussian) - CPU or GPU (if available)
% 平滑：可分离高斯卷积，支持 GPU（若可用）
% -----------------------------
rate_map_smooth = rate_map_raw; % default: no smoothing
if doSmooth
    bin_w = smooth_bin;
    smooth_sigma = - (bin_w^2) * 0.25 / log(0.5); % convert bin width to sigma heuristic
    if smooth_sigma > 0
        half = max(1, ceil(3 * smooth_sigma));    % kernel half-width (3 sigma rule)
        xv = -half:half;
        g1 = exp(-(double(xv).^2) / (2 * double(smooth_sigma)^2));
        g1 = g1 / sum(g1);                        % normalize kernel

        if useGPU && gpuDeviceCount > 0
            try
                rate_gpu = gpuArray(single(rate_map_raw)); % use single on GPU
                tmp = conv2(rate_gpu, reshape(single(g1),1,[]), 'same'); % convolve along columns
                rate_gpu = conv2(tmp, reshape(single(g1),[],1), 'same'); % then rows
                rate_map_smooth = gather(rate_gpu); % bring back to CPU
            catch ME
                warning('GPU smoothing failed (%s). Falling back to CPU. / GPU 平滑失败，回退到 CPU。', ME.message);
                tmp = conv2(rate_map_raw, reshape(g1,1,[]), 'same');
                rate_map_smooth = conv2(tmp, reshape(g1,[],1), 'same');
            end
        else
            % CPU separable convolution (fast and memory-friendly)
            tmp = conv2(rate_map_raw, reshape(g1,1,[]), 'same');
            rate_map_smooth = conv2(tmp, reshape(g1,[],1), 'same');
        end
    end
end

% keep NaN for info calculation but provide smoothed map for other uses
rate_map_for_output = rate_map_smooth;

% -----------------------------
% Spatial information (Skaggs bits/spike-like)
% 计算空间信息量（Skaggs 风格）
% -----------------------------
occ_vec = occupancy_time(:);
valid_occ = occ_vec > 0; % only bins that were actually visited / 仅考虑被访问的bin
if ~any(valid_occ)
    P = [];
else
    P = occ_vec(valid_occ) / sum(occ_vec(valid_occ)); % Pi occupancy probabilities
end

% overall mean activity rate R = total activity / total time visited
if isempty(P)
    R = 0;
else
    R = nansum(activity_map_raw(valid_occ)) / sum(occ_vec(valid_occ));
end

Ri = rate_map_smooth(:);
Ri = Ri(valid_occ); % align with Pi
idx_pos = ~isnan(Ri) & (Ri > 0); % bins with positive rate
if R == 0 || isempty(P) || ~any(idx_pos)
    info_obs = 0;
else
    info_obs = sum( P(idx_pos) .* (Ri(idx_pos) ./ R) .* log2( Ri(idx_pos) ./ R ) );
end

% -----------------------------
% Shuffle test (circular shift preserves temporal structure)
% 置换检验：spike 用连续时间 circular shift；CA 用帧级 circshift
% -----------------------------
nS = nShuffles;
info_shuf = zeros(nS, 1);

if nS > 0
    % precompute locals to avoid repeated capture overhead inside parfor
    valid_occ_local = valid_occ;
    P_local = P;
    R_local = R;
    sample_counts_local = sample_counts;
    occupancy_time_local = occupancy_time;
    x_edges_local = x_edges;
    y_edges_local = y_edges;
    nxny = [nx, ny];
    doSmooth_local = doSmooth;
    smooth_sigma_local = - (smooth_bin^2) * 0.25 / log(0.5);
    bin_pass_local = bin_pass;

    % CA precomputed bin-index map used in CA shuffles
    if data_type ~= "spike"
        inds = find(include_mask);
        xbin_full = nan(T,1);
        ybin_full = nan(T,1);
        xbin_full(inds) = xbin_idx;
        ybin_full(inds) = ybin_idx;
    end

    % create independent RNG seeds for reproducibility in parfor
    rng('shuffle');
    seeds = randi(intmax('int32'), nS, 1);

    if useParpool
        % Parallel loop - requires Parallel Computing Toolbox and a working pool
        parfor si = 1:nS
            % set per-iteration RNG stream to ensure independent randomness
            s = RandStream('mt19937ar', 'Seed', seeds(si));
            RandStream.setGlobalStream(s);

            if data_type == "spike"
                % continuous time circular shift for spike timestamps
                t_start = position_timelite(1);
                t_end = position_timelite(end);
                dur = t_end - t_start;
                shift = rand() * dur;
                spk_sh = mod(spike_times - t_start + shift, dur) + t_start;

                % map shifted spikes to positions using single interp + mask
                spk_x_all = interp1(position_timelite, positionX, spk_sh, 'linear', NaN);
                spk_y_all = interp1(position_timelite, positionY, spk_sh, 'linear', NaN);
                if ~isempty(intervalT0)
                    inMask = false(numel(spk_sh),1);
                    for kk = 1:numel(intervalT0)
                        inMask = inMask | (spk_sh >= intervalT0(kk) & spk_sh <= intervalT1(kk));
                    end
                else
                    inMask = false(numel(spk_sh),1);
                end
                spk_x_all(~inMask) = NaN;
                spk_y_all(~inMask) = NaN;

                xb = discretize(spk_x_all, x_edges_local);
                yb = discretize(spk_y_all, y_edges_local);
                valid_sp = ~isnan(xb) & ~isnan(yb);
                if any(valid_sp)
                    linIdx = sub2ind(nxny, xb(valid_sp), yb(valid_sp));
                    counts_lin = accumarray(linIdx, 1, [nxny(1) * nxny(2), 1]);
                    spk_map_sh = reshape(counts_lin, nxny);
                else
                    spk_map_sh = zeros(nxny);
                end
                activity_map_sh = spk_map_sh;
            else
                % CA mode shuffle: integer-frame circular shift preserves autocorrelation structure
                Tlen = numel(activity_vec);
                shift_frames = randi(Tlen) - 1;
                act_sh = circshift(activity_vec, shift_frames);
                valid_samples = ~isnan(xbin_full) & ~isnan(ybin_full) & ~isnan(act_sh);
                if any(valid_samples)
                    lin_idx = sub2ind(nxny, xbin_full(valid_samples), ybin_full(valid_samples));
                    sums = accumarray(lin_idx, double(act_sh(valid_samples)), [nxny(1) * nxny(2), 1]);
                    activity_map_sh = reshape(sums, nxny);
                else
                    activity_map_sh = zeros(nxny);
                end
            end

            % compute rate (activity/sec) for shuffle
            rate_sh = activity_map_sh ./ occupancy_time_local;
            rate_sh(~isfinite(rate_sh)) = NaN;
            if bin_pass_local > 0
                rate_sh(sample_counts_local < bin_pass_local) = 0;
            end

            % smoothing in parfor uses CPU separable conv for robustness
            if doSmooth_local && smooth_sigma_local > 0
                half = max(1, ceil(3 * smooth_sigma_local));
                xv = -half:half;
                g1 = exp(-(xv.^2) / (2 * smooth_sigma_local^2));
                g1 = g1 / sum(g1);
                tmp = conv2(rate_sh, reshape(g1,1,[]), 'same');
                rate_sh = conv2(tmp, reshape(g1,[],1), 'same');
            end

            % compute info for shuffled map
            Ri_sh = rate_sh(:);
            Ri_sh = Ri_sh(valid_occ_local);
            idx_pos_sh = ~isnan(Ri_sh) & (Ri_sh > 0);
            if R_local == 0 || ~any(idx_pos_sh)
                info_shuf(si) = 0;
            else
                info_shuf(si) = sum( P_local(idx_pos_sh) .* (Ri_sh(idx_pos_sh) ./ R_local) .* log2( Ri_sh(idx_pos_sh) ./ R_local ) );
            end
        end
    else
        % Serial loop (no parallel)
        for si = 1:nS
            s = RandStream('mt19937ar','Seed', seeds(si));
            RandStream.setGlobalStream(s);

            if data_type == "spike"
                t_start = position_timelite(1);
                t_end = position_timelite(end);
                dur = t_end - t_start;
                shift = rand() * dur;
                spk_sh = mod(spike_times - t_start + shift, dur) + t_start;
                spk_x_all = interp1(position_timelite, positionX, spk_sh, 'linear', NaN);
                spk_y_all = interp1(position_timelite, positionY, spk_sh, 'linear', NaN);
                if ~isempty(intervalT0)
                    inMask = false(numel(spk_sh),1);
                    for kk = 1:numel(intervalT0)
                        inMask = inMask | (spk_sh >= intervalT0(kk) & spk_sh <= intervalT1(kk));
                    end
                else
                    inMask = false(numel(spk_sh),1);
                end
                spk_x_all(~inMask) = NaN;
                spk_y_all(~inMask) = NaN;
                xb = discretize(spk_x_all, x_edges_local);
                yb = discretize(spk_y_all, y_edges_local);
                valid_sp = ~isnan(xb) & ~isnan(yb);
                if any(valid_sp)
                    linIdx = sub2ind(nxny, xb(valid_sp), yb(valid_sp));
                    counts_lin = accumarray(linIdx, 1, [nxny(1) * nxny(2), 1]);
                    spk_map_sh = reshape(counts_lin, nxny);
                else
                    spk_map_sh = zeros(nxny);
                end
                activity_map_sh = spk_map_sh;
            else
                Tlen = numel(activity_vec);
                shift_frames = randi(Tlen) - 1;
                act_sh = circshift(activity_vec, shift_frames);
                valid_samples = ~isnan(xbin_full) & ~isnan(ybin_full) & ~isnan(act_sh);
                if any(valid_samples)
                    lin_idx = sub2ind(nxny, xbin_full(valid_samples), ybin_full(valid_samples));
                    sums = accumarray(lin_idx, double(act_sh(valid_samples)), [nxny(1) * nxny(2), 1]);
                    activity_map_sh = reshape(sums, nxny);
                else
                    activity_map_sh = zeros(nxny);
                end
            end

            rate_sh = activity_map_sh ./ occupancy_time_local;
            rate_sh(~isfinite(rate_sh)) = NaN;
            if bin_pass_local > 0
                rate_sh(sample_counts_local < bin_pass_local) = 0;
            end
            if doSmooth_local && smooth_sigma_local > 0
                half = max(1, ceil(3 * smooth_sigma_local));
                xv = -half:half;
                g1 = exp(-(xv.^2) / (2 * smooth_sigma_local^2));
                g1 = g1 / sum(g1);
                tmp = conv2(rate_sh, reshape(g1,1,[]), 'same');
                rate_sh = conv2(tmp, reshape(g1,[],1), 'same');
            end
            Ri_sh = rate_sh(:);
            Ri_sh = Ri_sh(valid_occ_local);
            idx_pos_sh = ~isnan(Ri_sh) & (Ri_sh > 0);
            if R_local == 0 || ~any(idx_pos_sh)
                info_shuf(si) = 0;
            else
                info_shuf(si) = sum( P_local(idx_pos_sh) .* (Ri_sh(idx_pos_sh) ./ R_local) .* log2( Ri_sh(idx_pos_sh) ./ R_local ) );
            end
        end
    end
end

% compute p-value (robust, +1 / +1)
if nS > 0
    pval = (sum(info_shuf >= info_obs) + 1) / (nS + 1);
else
    pval = NaN;
end

% -----------------------------
% Detect place fields on smoothed map
% 在平滑后的地图上检测 place field
% -----------------------------
map_for_field = rate_map_smooth;
peakRate = max(map_for_field(~isnan(map_for_field)), [], 'all');
if isempty(peakRate) || isnan(peakRate)
    peakRate = 0;
end
thresh_val = thr_frac * peakRate; % threshold value to define candidate field mask
bw_mask = (map_for_field >= thresh_val) & ~isnan(map_for_field);

% Connected components (8-connectivity)
CC = bwconncomp(bw_mask, 8);
if CC.NumObjects > 0
    props = regionprops(CC, map_for_field, 'Area', 'PixelIdxList', 'MaxIntensity', 'Centroid');
else
    props = [];
end

% build fields struct with peakRate, area_bins, centroid_x/y, mask
fields = struct('peakRate', {}, 'area_bins', {}, 'centroid_x', {}, 'centroid_y', {}, 'mask', {});
x_centers = x_edges(1:end-1) + diff(x_edges)/2;
y_centers = y_edges(1:end-1) + diff(y_edges)/2;
for k = 1:CC.NumObjects
    pix = CC.PixelIdxList{k};
    area_bins = numel(pix);
    peak_r = props(k).MaxIntensity;
    [ix, iy] = ind2sub(size(map_for_field), pix); % ix->x-bin idx, iy->y-bin idx (matrix coords)
    cx = mean(x_centers(ix)); % centroid x in original coordinates
    cy = mean(y_centers(iy)); % centroid y in original coordinates
    mask = false(size(map_for_field));
    mask(pix) = true;
    fields(k).peakRate = peak_r;
    fields(k).area_bins = area_bins;
    fields(k).centroid_x = cx;
    fields(k).centroid_y = cy;
    fields(k).mask = mask;
end

% filter fields by minimal area and peak rate
valid_fields = [];
for k = 1:numel(fields)
    if fields(k).area_bins >= min_field_bins && fields(k).peakRate >= peak_rate_min
        valid_fields = [valid_fields, fields(k)]; %#ok<AGROW>
    end
end

% final place cell decision: info threshold, p-value threshold, and at least one valid field
is_place = (info_obs >= info_thresh) && (pval <= p_thresh) && (~isempty(valid_fields));

% -----------------------------
% Populate results struct (same fields as original + parameter flags)
% 填充输出结构（兼容原输出并增加参数标志）
% -----------------------------
results = struct();
results.is_place_cell = logical(is_place);
results.info = info_obs;
results.info_shuffles = info_shuf;
results.p_value = pval;
results.rate_map_raw = rate_map_raw;
results.rate_map_smooth = rate_map_smooth;
results.spike_map = spk_map;                % NaN for CA mode
results.activity_map_raw = activity_map_raw; % counts or summed activity
results.occupancy_time = occupancy_time;
results.sample_counts = sample_counts;
results.x_edges = x_edges;
results.y_edges = y_edges;
results.fields_all = fields;
results.fields_valid = valid_fields;
results.peakRate = peakRate;
params = p.Results;
params.bin_size = bin_size;
params.useParpool = useParpool;
params.useGPU = useGPU;
params.useSingle = useSingle;
results.parameters = params;

% Verbose printout for quick summary
if verbose
    fprintf('Mode: %s. Spatial info: %.3f bits, p=%.4f (nShuffles=%d), detected %d valid fields. (useParpool=%d,useGPU=%d,useSingle=%d)\n', ...
        char(data_type), info_obs, pval, nS, numel(valid_fields), double(useParpool), double(useGPU), double(useSingle));
    % 打印简短摘要（中英）:
    fprintf('模式: %s。空间信息: %.3f bits, p=%.4f (置换次数=%d)，检测到 %d 个有效场。（useParpool=%d,useGPU=%d,useSingle=%d）\n', ...
        char(data_type), info_obs, pval, nS, numel(valid_fields), double(useParpool), double(useGPU), double(useSingle));
end

end
