function [H, x_edges, y_edges, H_intervals] = occupancy(x, y, varargin)
% occupancy_hist2d_interval_v5
% 计算轨迹二维占用直方图（timeinterval 为 T×n，每列应为单一连续块）
%
% USAGE:
%   [H, xe, ye] = occupancy_hist2d_interval_v5(x, y)
%   [H, xe, ye] = occupancy_hist2d_interval_v5(x, y, 'timeinterval', TI, 'fps',30)
%   [H, xe, ye, Hints] = occupancy_hist2d_interval_v5(..., 'per_interval', true)
%
% NAME-VALUE:
%   'bin_size'     (default 10)
%   'fps'          (default 30)
%   'timeinterval' - T×n logical matrix (preferred), or cell/numeric as before
%                    **Each column is expected to have one contiguous block of 1's.**
%   'speed'        - speed vector (T×1) if using 't_threshold'
%   't_threshold'  - scalar threshold for speed (exclude speed < threshold)
%   'normalize'    - logical (default false) normalize H to sum(H)=1
%   'per_interval' - logical (default false), if true also return H_intervals (rows: y bins, cols: x bins, 3rd dim: interval index)
%
% OUTPUT:
%   H            - occupancy map (seconds), union over all given intervals (or all time if no TI)
%   x_edges,y_edges
%   H_intervals  - optional (if per_interval true): H_intervals(:,:,k) is occupancy (s) for k-th interval

% parse inputs
p = inputParser;
addParameter(p,'bin_size',10,@(v)isnumeric(v)&&isscalar(v)&&v>0);
addParameter(p,'fps',30,@(v)isnumeric(v)&&isscalar(v)&&v>0);
addParameter(p,'timeinterval',[], @(v) (islogical(v)||iscell(v)||isnumeric(v)||isempty(v)));
addParameter(p,'speed',[], @(v)isnumeric(v)&&(isvector(v)||isempty(v)));
addParameter(p,'t_threshold',[], @(v) isempty(v) || (isnumeric(v)&&isscalar(v)));
addParameter(p,'normalize',false, @(v)islogical(v)||(isnumeric(v)&&isscalar(v)));
addParameter(p,'per_interval',false, @(v) islogical(v) || (isnumeric(v)&&isscalar(v)));
parse(p,varargin{:});
bin_size   = p.Results.bin_size;
fps        = p.Results.fps;
TI_in      = p.Results.timeinterval;
spd        = p.Results.speed;
t_threshold= p.Results.t_threshold;
doNorm     = logical(p.Results.normalize);
perInterval= logical(p.Results.per_interval);

% ensure column vectors
x = x(:); y = y(:);
T = numel(x);
if numel(y) ~= T, error('x and y must have same length (T).'); end

% initial include mask (exclude NaNs)
includeMask = ~(isnan(x) | isnan(y));

% process TI: prefer logical T x n
TI_proc = []; % will be logical T x n or []
if ~isempty(TI_in)
    if islogical(TI_in)
        if size(TI_in,1) ~= T
            error('timeinterval must have T rows equal to length(x).');
        end
        TI_proc = TI_in;
    elseif iscell(TI_in)
        % convert cell of index vectors to logical matrix
        ncols = numel(TI_in);
        TI_proc = false(T, ncols);
        for k=1:ncols
            idx = TI_in{k}(:);
            idx = idx(idx>=1 & idx<=T);
            TI_proc(idx,k) = true;
        end
    elseif isnumeric(TI_in)
        % numeric matrix: non-zero indicates inclusion; expect T x n (or transpose)
        if size(TI_in,1) ~= T
            if size(TI_in,2) == T
                TI_in = TI_in';
            else
                error('timeinterval numeric matrix must have T rows equal to length(x).');
            end
        end
        TI_proc = TI_in ~= 0;
    else
        error('Unsupported timeinterval type.');
    end
    % Validate columns are (or pick longest contiguous block)
    n_intervals = size(TI_proc,2);
    for k = 1:n_intervals
        col = TI_proc(:,k);
        if ~any(col)
            warning('timeinterval column %d is empty (all false). It will be ignored.', k);
            continue;
        end
        % find contiguous true segments
        d = diff([0; col; 0]);
        starts = find(d==1);
        ends   = find(d==-1)-1;
        nseg = numel(starts);
        if nseg > 1
            % multiple segments: warn and keep the longest
            seglens = ends - starts + 1;
            [~, idxmax] = max(seglens);
            % build new column that keeps only longest segment
            newcol = false(T,1);
            newcol(starts(idxmax):ends(idxmax)) = true;
            TI_proc(:,k) = newcol;
            warning('timeinterval column %d has %d disjoint segments: keeping longest segment [%d:%d].', k, nseg, starts(idxmax), ends(idxmax));
        end
    end
    % union mask
    includeMask = includeMask & any(TI_proc,2);
end

% apply speed threshold if provided
if ~isempty(t_threshold)
    if isempty(spd), error('When specifying t_threshold, you must supply speed vector.'); end
    spd = spd(:);
    if numel(spd) ~= T, error('speed must have length T.'); end
    includeMask = includeMask & (spd >= t_threshold);
end

% if no points remain
if ~any(includeMask)
    warning('No valid points after applying NaN/timeinterval/speed filters. Returning empty histogram.');
    x_edges = [0 bin_size];
    y_edges = [0 bin_size];
    H = zeros(length(y_edges)-1, length(x_edges)-1);
    if perInterval, H_intervals = []; end
    return;
end

% compute edges from included points
x_use = x(includeMask);
y_use = y(includeMask);
x_min = min(x_use); x_max = max(x_use);
y_min = min(y_use); y_max = max(y_use);
if x_min == x_max
    x_min = x_min - bin_size/2; x_max = x_max + bin_size/2;
end
if y_min == y_max
    y_min = y_min - bin_size/2; y_max = y_max + bin_size/2;
end
x_edges = x_min:bin_size:(x_max+eps);
y_edges = y_min:bin_size:(y_max+eps);

% compute union histogram (counts -> seconds)
H_counts = histcounts2( x_use,y_use, x_edges, y_edges);
H = H_counts / fps;

% optionally compute per-interval histograms
H_intervals = [];
if perInterval && ~isempty(TI_proc)
    nI = size(TI_proc,2);
    H_intervals = zeros(size(H,1), size(H,2), nI);
    for k = 1:nI
        idxk = TI_proc(:,k);
        % ensure also exclude NaNs and speed-thresholded points
        idxk = idxk & ~(isnan(x) | isnan(y));
        if ~isempty(t_threshold)
            idxk = idxk & (spd >= t_threshold);
        end
        if ~any(idxk)
            Hk = zeros(size(H));
        else
            Hk_counts = histcounts2(x(idxk),y(idxk), x_edges,  y_edges);
            Hk = Hk_counts / fps;
        end
        H_intervals(:,:,k) = Hk;
    end
end

% normalize if requested (after conversion to time)
if doNorm
    s = sum(H(:));
    if s>0, H = H / s; end
    if perInterval && ~isempty(H_intervals)
        for k=1:size(H_intervals,3)
            sk = sum(H_intervals(:,:,k),'all');
            if sk>0, H_intervals(:,:,k) = H_intervals(:,:,k) / sk; end
        end
    end
end
end
