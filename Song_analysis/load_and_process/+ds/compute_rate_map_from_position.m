function [rate_map, x_edges, y_edges] = compute_rate_map_from_position(...
    positionX, positionY, position_timelite, spike_times, inIntervals, bin_size, varargin)
% COMPUTE_RATE_MAP_FROM_POSITION
%   Compute 2D firing rate map from spike times & positions, optionally smoothed.
%
% USAGE:
%   rate_map = compute_rate_map_from_position(posX,posY,posT,spike_times,inIntervals,bin_size)
%   rate_map = compute_rate_map_from_position(..., 'smooth', true, 'smooth_bin', 2, 'bin_pass', 5)
%
% REQUIRED:
%   positionX, positionY, position_timelite (Tx1)
%   spike_times (Sx1)
%   inIntervals (T x n logical matrix or [])
%   bin_size (scalar)
%
% NAME-VALUE (defaults):
%   'smooth' (false)            - whether to Gaussian-smooth the rate map
%   'smooth_bin' (2)            - smoothing scale in bins (only used if smooth==true)
%   'map_nan_to_zero' (true)    - convert NaN (no occupancy) to 0 in final output
%   'bin_pass' (0)              - minimum sample-count in a bin to consider its rate valid.
%                                 If sample_count < bin_pass, set rate to 0.
%
% OUTPUT:
%   rate_map (nx x ny)          - smoothed or raw rate map per 'smooth' setting
%   x_edges, y_edges            - bin edges used

%% parse inputs
p = inputParser;
addParameter(p,'smooth',false,@(v)islogical(v)||isnumeric(v));
addParameter(p,'smooth_bin',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'map_nan_to_zero',true,@(v)islogical(v)||isnumeric(v)); % default true
addParameter(p,'bin_pass',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});
doSmooth = logical(p.Results.smooth);
smooth_bin = p.Results.smooth_bin;
map_nan_to_zero = logical(p.Results.map_nan_to_zero);
bin_pass = p.Results.bin_pass;

%% basic checks
positionX = positionX(:);
positionY = positionY(:);
position_timelite = position_timelite(:);
spike_times = spike_times(:);
T = numel(positionX);
if numel(positionY) ~= T || numel(position_timelite) ~= T
    error('positionX, positionY, and position_timelite must have same length.');
end
if nargin < 5 || isempty(inIntervals)
    inIntervals = true(T,1);
end
if ~islogical(inIntervals)
    inIntervals = inIntervals ~= 0;
end
if size(inIntervals,1) ~= T
    error('inIntervals must have T rows.');
end
if nargin < 6 || isempty(bin_size)
    error('You must specify bin_size.');
end

%% sample durations dt (seconds)
dt = [diff(position_timelite); median(diff(position_timelite))];
if ~any(dt>0)
    dt(:) = 1;
end

%% define bin edges from valid positions within intervals
valid_pos = ~isnan(positionX) & ~isnan(positionY);
mask_any = any(inIntervals,2);
mask_edges = valid_pos & mask_any;
if ~any(mask_edges); mask_edges = valid_pos; end

xmin = min(positionX(mask_edges)); xmax = max(positionX(mask_edges));
ymin = min(positionY(mask_edges)); ymax = max(positionY(mask_edges));
pad = bin_size/2;
x_edges = xmin-pad : bin_size : xmax+pad;
y_edges = ymin-pad : bin_size : ymax+pad;
nx = numel(x_edges)-1;
ny = numel(y_edges)-1;

%% occupancy: compute sample counts and total occupancy time per bin
% include_mask = samples counted for occupancy (union of intervals & valid positions)
include_mask = valid_pos & mask_any;
% sample counts per bin (integer number of position samples)
if any(include_mask)
    sample_counts = histcounts2(positionX(include_mask), positionY(include_mask), x_edges, y_edges);
else
    sample_counts = zeros(nx, ny);
end
% occupancy_time (seconds) computed by summing dt for each included sample
% We'll compute with histcounts2 over weighted values: since histcounts2 doesn't accept weights,
% accumulate via bin indices for efficiency:
occupancy_time = zeros(nx, ny);
if any(include_mask)
    [~, ~, xbin_idx] = histcounts(positionX(include_mask), x_edges);
    [~, ~, ybin_idx] = histcounts(positionY(include_mask), y_edges);
    inds = find(include_mask);
    % iterate over included samples (vectorized alternative could use accumarray)
    for k = 1:numel(inds)
        i = inds(k);
        xb = xbin_idx(k); yb = ybin_idx(k);
        if xb>=1 && xb<=nx && yb>=1 && yb<=ny
            occupancy_time(xb,yb) = occupancy_time(xb,yb) + dt(i);
        end
    end
end

%% Map spikes -> positions (per contiguous segment per interval)
S = numel(spike_times);
spike_x = NaN(S,1); spike_y = NaN(S,1);
[spk_sort, sortIdx] = sort(spike_times);
nCols = size(inIntervals,2);
for c = 1:nCols
    col = inIntervals(:,c);
    if ~any(col), continue; end
    d = diff([0; col(:); 0]);
    starts = find(d==1); ends = find(d==-1)-1;
    for si = 1:numel(starts)
        sidx = starts(si); eidx = ends(si);
        t_seg = position_timelite(sidx:eidx);
        x_seg = positionX(sidx:eidx);
        y_seg = positionY(sidx:eidx);
        % fill internal NaNs by linear interpolation when possible
        if sum(~isnan(x_seg)) >= 2
            x_seg = interp1(find(~isnan(x_seg)), x_seg(~isnan(x_seg)), (1:numel(x_seg))', 'linear', NaN);
        end
        if sum(~isnan(y_seg)) >= 2
            y_seg = interp1(find(~isnan(y_seg)), y_seg(~isnan(y_seg)), (1:numel(y_seg))', 'linear', NaN);
        end
        % find spikes in segment
        t0 = t_seg(1); t1 = t_seg(end);
        idxs_spk_sorted = find(spk_sort >= t0 & spk_sort <= t1);
        if isempty(idxs_spk_sorted), continue; end
        t_in = spk_sort(idxs_spk_sorted);
        x_in = interp1(t_seg, x_seg, t_in, 'linear', NaN);
        y_in = interp1(t_seg, y_seg, t_in, 'linear', NaN);
        valid_sp = ~isnan(x_in) & ~isnan(y_in);
        if any(valid_sp)
            orig_idx = sortIdx(idxs_spk_sorted(valid_sp));
            spike_x(orig_idx) = x_in(valid_sp);
            spike_y(orig_idx) = y_in(valid_sp);
        end
    end
end

%% spike_map (counts per bin)
valid_sp = ~isnan(spike_x) & ~isnan(spike_y);
if any(valid_sp)
    spk_map = histcounts2(spike_x(valid_sp), spike_y(valid_sp), x_edges, y_edges);
else
    spk_map = zeros(nx, ny);
end

%% raw rate map (spikes per second), NaN where occupancy_time == 0
rate_map = spk_map ./ occupancy_time;
rate_map(~isfinite(rate_map)) = 0;

%% apply bin_pass threshold (sample count threshold)
if bin_pass > 0
    low_sample_mask = sample_counts < bin_pass; % logical matrix
    % set those bins to 0 as requested
    rate_map(low_sample_mask) = 0;
end

%% optional smoothing (applied after bin_pass thresholding)
if doSmooth
    bin_w = smooth_bin;
    smooth_sigma = - (bin_w^2) * 0.25 / log(0.5);
    if smooth_sigma <= 0
        % nothing
    else
        try
            rate_map = imgaussfilt(rate_map, smooth_sigma, 'FilterSize', max(3,2*ceil(3*smooth_sigma)+1));
        catch
            % fallback conv2 gaussian
            sz = max(3,2*ceil(3*smooth_sigma)+1);
            half = floor(sz/2);
            xv = -half:half;
            g1 = exp(-(xv.^2)/(2*smooth_sigma^2)); g1 = g1 / sum(g1);
            G = g1' * g1;
            pad = floor(sz/2);
            tmp = padarray(rate_map, [pad pad], 'symmetric');
            tmp2 = conv2(tmp, G, 'same');
            rate_map = tmp2(pad+1:end-pad, pad+1:end-pad);
        end
    end
end

%% convert NaN -> 0 if requested (default true)
if map_nan_to_zero
    rate_map(isnan(rate_map)) = 0;
end

end
