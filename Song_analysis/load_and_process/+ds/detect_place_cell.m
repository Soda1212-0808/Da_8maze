function results = detect_place_cell(positionX, positionY, position_timelite, spike_times, inIntervals, bin_size, varargin)
% DETECT_PLACE_CELL  detect place fields and test whether a neuron is a place cell.
%
% USAGE:
%   results = detect_place_cell(positionX, positionY, position_timelite, spike_times, inIntervals, bin_size, ...)
%
% REQUIRED:
%   positionX, positionY      - Tx1
%   position_timelite         - Tx1 timestamps (seconds)
%   spike_times               - Sx1 (seconds)
%   inIntervals               - T x n logical matrix (or [])
%   bin_size                  - scalar (same units as positionX/Y)
%
% NAME-VALUE OPTIONS (sensible defaults):
%   'smooth' (false)              - whether to Gaussian-smooth the rate map
%   'smooth_bin' (2)              - smoothing scale in bins (used when smooth==true)
%   'bin_pass' (0)                - minimum sample-count per bin to consider valid
%   'min_t' ([])                  - if provided, use min_t (s) to compute bin_pass as ceil(min_t * frame_rate)
%   'frame_rate' ([])             - if min_t provided, required to compute bin_pass (Hz)
%   'nShuffles' (1000)            - number of circular-shift shuffles for significance
%   'info_thresh' (0.5)           - bits/spike threshold (common heuristic)
%   'p_thresh' (0.05)             - significance threshold for shuffle p-value
%   'min_field_bins' (3)          - minimal number of connected bins for a valid field
%   'peak_rate_min' (1)           - minimal peak firing rate (Hz) required for a valid field
%   'field_threshold_fraction' (0.2) - threshold fraction of peak to define field (default 0.2)
%   'verbose' (true)
%
% OUTPUT (structure 'results'):
%   results.is_place_cell (logical)
%   results.info              - observed spatial information (bits/spike)
%   results.p_value           - p-value from shuffling
%   results.rate_map_raw      - raw rate map (nx x ny)
%   results.rate_map_smooth   - smoothed rate map (or same as raw if smooth=false)
%   results.spike_map         - spike counts per bin
%   results.occupancy_time    - seconds per bin
%   results.x_edges, y_edges
%   results.fields            - array of structs (peakRate, area_bins, centroid_x, centroid_y, mask)
%   results.parameters        - parameters used
%
% NOTE: matrix orientation: first dimension = x bins, second = y bins (consistent with histcounts2 usage).

% parse name-value
p = inputParser;
addParameter(p,'smooth',false,@(v)islogical(v)||isnumeric(v));
addParameter(p,'smooth_bin',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'bin_pass',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'min_t',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(p,'frame_rate',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'nShuffles',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'info_thresh',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'p_thresh',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'min_field_bins',3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'peak_rate_min',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'field_threshold_fraction',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'verbose',true,@(v)islogical(v)||isnumeric(v));
parse(p,varargin{:});
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

% auto compute bin_pass from min_t if requested
if ~isempty(min_t)
    if isempty(frame_rate)
        error('frame_rate must be provided when min_t is set.');
    end
    bin_pass = ceil(min_t * frame_rate);
end

% ---- Basic checks / reshape ----
positionX = positionX(:);
positionY = positionY(:);
position_timelite = position_timelite(:);
spike_times = spike_times(:);
T = numel(positionX);
if numel(positionY) ~= T || numel(position_timelite) ~= T
    error('positionX/Y/timelite must have same length T.');
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

% ---- compute occupancy and spike_map (same logic as compute_rate_map_from_position) ----
% build edges
valid_pos = ~isnan(positionX) & ~isnan(positionY);
mask_any = any(inIntervals,2);
mask_edges = valid_pos & mask_any;
if ~any(mask_edges), mask_edges = valid_pos; end

xmin = min(positionX(mask_edges)); xmax = max(positionX(mask_edges));
ymin = min(positionY(mask_edges)); ymax = max(positionY(mask_edges));
pad = bin_size/2;
x_edges = xmin-pad : bin_size : xmax+pad;
y_edges = ymin-pad : bin_size : ymax+pad;
nx = numel(x_edges)-1;
ny = numel(y_edges)-1;

% dt
dt = [diff(position_timelite); median(diff(position_timelite))];
if ~any(dt>0)
    dt(:)=1;
end

% occupancy sample counts (and occupancy_time)
include_mask = valid_pos & mask_any;
if any(include_mask)
    sample_counts = histcounts2(positionX(include_mask), positionY(include_mask), x_edges, y_edges);
else
    sample_counts = zeros(nx, ny);
end
occupancy_time = zeros(nx, ny);
if any(include_mask)
    [~,~,xbin_idx] = histcounts(positionX(include_mask), x_edges);
    [~,~,ybin_idx] = histcounts(positionY(include_mask), y_edges);
    inds = find(include_mask);
    for k = 1:numel(inds)
        i = inds(k);
        xb = xbin_idx(k); yb = ybin_idx(k);
        if xb>=1 && xb<=nx && yb>=1 && yb<=ny
            occupancy_time(xb,yb) = occupancy_time(xb,yb) + dt(i);
        end
    end
end

% map spikes to positions (same per-interval segmentation as before)
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
        if sum(~isnan(x_seg)) >= 2
            x_seg = interp1(find(~isnan(x_seg)), x_seg(~isnan(x_seg)), (1:numel(x_seg))','linear',NaN);
        end
        if sum(~isnan(y_seg)) >= 2
            y_seg = interp1(find(~isnan(y_seg)), y_seg(~isnan(y_seg)), (1:numel(y_seg))','linear',NaN);
        end
        t0 = t_seg(1); t1 = t_seg(end);
        idxs_spk_sorted = find(spk_sort >= t0 & spk_sort <= t1);
        if isempty(idxs_spk_sorted), continue; end
        t_in = spk_sort(idxs_spk_sorted);
        x_in = interp1(t_seg, x_seg, t_in, 'linear', NaN);
        y_in = interp1(t_seg, y_seg, t_in, 'linear', NaN);
        v = ~isnan(x_in) & ~isnan(y_in);
        if any(v)
            orig_idx = sortIdx(idxs_spk_sorted(v));
            spike_x(orig_idx) = x_in(v);
            spike_y(orig_idx) = y_in(v);
        end
    end
end

% spike_map
valid_sp = ~isnan(spike_x) & ~isnan(spike_y);
if any(valid_sp)
    spk_map = histcounts2(spike_x(valid_sp), spike_y(valid_sp), x_edges, y_edges);
else
    spk_map = zeros(nx, ny);
end

% raw rate map
rate_map_raw = spk_map ./ occupancy_time;
rate_map_raw(~isfinite(rate_map_raw)) = NaN;

% apply bin_pass (sample-count threshold)
if bin_pass > 0
    low_mask = sample_counts < bin_pass;
    rate_map_raw(low_mask) = 0;  % as you requested: low-visit bins considered 0
end

% smoothing if requested
rate_map_smooth = rate_map_raw;
if doSmooth
    bin_w = smooth_bin;
    smooth_sigma = - (bin_w^2)*0.25 / log(0.5);
    if smooth_sigma > 0
        try
            rate_map_smooth = imgaussfilt(rate_map_raw, smooth_sigma, 'FilterSize', max(3,2*ceil(3*smooth_sigma)+1));
        catch
            sz = max(3,2*ceil(3*smooth_sigma)+1);
            half = floor(sz/2);
            xv = -half:half;
            g1 = exp(-(xv.^2)/(2*smooth_sigma^2)); g1 = g1/sum(g1);
            G = g1'*g1;
            pad = floor(sz/2);
            tmp = padarray(rate_map_raw,[pad pad],'symmetric');
            tmp2 = conv2(tmp,G,'same');
            rate_map_smooth = tmp2(pad+1:end-pad,pad+1:end-pad);
        end
    end
end

% convert NaN->0 for any visualization or downstream steps? we keep NaN for info calc
rate_map_for_output = rate_map_smooth;

% compute spatial information (Skaggs bits/spike)
% use occupancy probabilities Pi from occupancy_time (exclude zero occupancy bins)
occ_vec = occupancy_time(:);
valid_occ = occ_vec > 0;
P = occ_vec(valid_occ) / sum(occ_vec(valid_occ));            % Pi
R = nansum(spk_map(:)) / sum(occ_vec(valid_occ));            % overall mean firing rate (spikes/sec)
Ri = rate_map_smooth(:);
Ri = Ri(valid_occ);  % align with Pi
% avoid negative or zero issues: define only bins where Ri>0
idx_pos = ~isnan(Ri) & (Ri > 0);
if R == 0 || sum(P)<=0 || ~any(idx_pos)
    info_obs = 0;
else
    info_obs = sum( P(idx_pos) .* (Ri(idx_pos)./R) .* log2( Ri(idx_pos)./R ) );
end

% shuffle test: circular time-shift of spike_times (preserve ISI structure)
nS = nShuffles;
info_shuf = zeros(nS,1);
if nS > 0
    t_start = position_timelite(1);
    t_end = position_timelite(end);
    dur = t_end - t_start;
    rng('shuffle');
    parfor_or_for = @forloop; % helper to allow optional parfor later
    % We'll use for (not parfor) to avoid requiring Parallel Toolbox; but user can modify
    for si = 1:nS
        shift = rand()*dur; % shift in seconds
        spk_sh = mod(spike_times - t_start + shift, dur) + t_start;
        % remap spikes to rate map using same procedure but only need info -> compute spike_map_sh
        % For speed we map spike_sh to bins by linear interpolation on position_timelite segments as above
        % Reuse same mapping logic but inline simplified: use interp over entire concatenated intervals? we use same per-interval mapping
        spike_x_sh = NaN(S,1); spike_y_sh = NaN(S,1);
        [spk_sort_sh, sortIdx_sh] = sort(spk_sh);
        for c = 1:nCols
            col = inIntervals(:,c);
            if ~any(col), continue; end
            d = diff([0; col(:); 0]);
            starts = find(d==1); ends = find(d==-1)-1;
            for ssi = 1:numel(starts)
                si1=starts(ssi); si2=ends(ssi);
                t_seg = position_timelite(si1:si2);
                x_seg = positionX(si1:si2);
                y_seg = positionY(si1:si2);
                if sum(~isnan(x_seg))>=2
                    x_seg = interp1(find(~isnan(x_seg)), x_seg(~isnan(x_seg)), (1:numel(x_seg))','linear',NaN);
                end
                if sum(~isnan(y_seg))>=2
                    y_seg = interp1(find(~isnan(y_seg)), y_seg(~isnan(y_seg)), (1:numel(y_seg))','linear',NaN);
                end
                t0 = t_seg(1); t1 = t_seg(end);
                idxs_sorted = find(spk_sort_sh >= t0 & spk_sort_sh <= t1);
                if isempty(idxs_sorted), continue; end
                t_in = spk_sort_sh(idxs_sorted);
                x_in = interp1(t_seg, x_seg, t_in, 'linear', NaN);
                y_in = interp1(t_seg, y_seg, t_in, 'linear', NaN);
                v = ~isnan(x_in)&~isnan(y_in);
                if any(v)
                    orig_idx2 = sortIdx_sh(idxs_sorted(v));
                    spike_x_sh(orig_idx2) = x_in(v);
                    spike_y_sh(orig_idx2) = y_in(v);
                end
            end
        end
        % spike_map_sh
        valid_sh = ~isnan(spike_x_sh) & ~isnan(spike_y_sh);
        if any(valid_sh)
            spk_map_sh = histcounts2(spike_x_sh(valid_sh), spike_y_sh(valid_sh), x_edges, y_edges);
        else
            spk_map_sh = zeros(nx,ny);
        end
        % if bin_pass applied earlier, we apply same threshold to shuf maps too
        rate_sh = spk_map_sh ./ occupancy_time;
        rate_sh(~isfinite(rate_sh)) = NaN;
        if bin_pass > 0
            rate_sh(sample_counts < bin_pass) = 0;
        end
        if doSmooth
            % smooth same way
            bin_w = smooth_bin;
            smooth_sigma = - (bin_w^2)*0.25 / log(0.5);
            if smooth_sigma > 0
                try
                    rate_sh = imgaussfilt(rate_sh, smooth_sigma, 'FilterSize', max(3,2*ceil(3*smooth_sigma)+1));
                catch
                    sz = max(3,2*ceil(3*smooth_sigma)+1);
                    half = floor(sz/2);
                    xv = -half:half;
                    g1 = exp(-(xv.^2)/(2*smooth_sigma^2)); g1=g1/sum(g1);
                    G = g1'*g1;
                    pad = floor(sz/2);
                    tmp = padarray(rate_sh,[pad pad],'symmetric');
                    tmp2 = conv2(tmp,G,'same');
                    rate_sh = tmp2(pad+1:end-pad,pad+1:end-pad);
                end
            end
        end
        % compute info for shuffle
        Ri_sh = rate_sh(:);
        Ri_sh = Ri_sh(valid_occ);
        idx_pos_sh = ~isnan(Ri_sh) & (Ri_sh > 0);
        if R == 0 || ~any(idx_pos_sh)
            info_shuf(si) = 0;
        else
            info_shuf(si) = sum( P(idx_pos_sh) .* (Ri_sh(idx_pos_sh)./R) .* log2( Ri_sh(idx_pos_sh)./R ) );
        end
    end
else
    info_shuf = zeros(0,1);
end

% compute p-value
if nS > 0
    pval = (sum(info_shuf >= info_obs) + 1) / (nS + 1); % +1 for robust estimate
else
    pval = NaN;
end

% detect fields on the smoothed map (or raw if not smoothed)
map_for_field = rate_map_smooth;
% define threshold = thr_frac * peakRate
peakRate = max(map_for_field(~isnan(map_for_field)), [], 'all');
if isempty(peakRate) || isnan(peakRate)
    peakRate = 0;
end
thresh_val = thr_frac * peakRate;
bw_mask = (map_for_field >= thresh_val) & ~isnan(map_for_field);

CC = bwconncomp(bw_mask, 8); % 8-connectivity
props = regionprops(CC, map_for_field, 'Area', 'PixelIdxList', 'MaxIntensity', 'Centroid');

fields = struct('peakRate',{},'area_bins',{},'centroid_x',{},'centroid_y',{},'mask',{});
% compute bin centers
x_centers = x_edges(1:end-1) + diff(x_edges)/2;
y_centers = y_edges(1:end-1) + diff(y_edges)/2;
for k = 1:CC.NumObjects
    pix = CC.PixelIdxList{k}; % linear indices into map (column-major)
    area_bins = numel(pix);
    peak_r = props(k).MaxIntensity;
    % centroid in matrix coordinates -> convert to x,y via sub2ind
    [ix, iy] = ind2sub(size(map_for_field), pix); % ix->x-bin idx, iy->y-bin idx
    cx = mean(x_centers(ix));
    cy = mean(y_centers(iy));
    mask = false(size(map_for_field));
    mask(pix) = true;
    fields(k).peakRate = peak_r;
    fields(k).area_bins = area_bins;
    fields(k).centroid_x = cx;
    fields(k).centroid_y = cy;
    fields(k).mask = mask;
end

% filter fields by min bins and peak rate threshold
valid_fields = [];
for k=1:numel(fields)
    if fields(k).area_bins >= min_field_bins && fields(k).peakRate >= peak_rate_min
        valid_fields = [valid_fields, fields(k)];
    end
end

% final place cell decision
is_place = (info_obs >= info_thresh) && (pval <= p_thresh) && (~isempty(valid_fields));

% fill results
results = struct();
results.is_place_cell = logical(is_place);
results.info = info_obs;
results.info_shuffles = info_shuf;
results.p_value = pval;
results.rate_map_raw = rate_map_raw;
results.rate_map_smooth = rate_map_smooth;
results.spike_map = spk_map;
results.occupancy_time = occupancy_time;
results.sample_counts = sample_counts;
results.x_edges = x_edges;
results.y_edges = y_edges;
results.fields_all = fields;
results.fields_valid = valid_fields;
results.peakRate = peakRate;
results.parameters = p.Results;
results.parameters.bin_size = bin_size;

if verbose
    fprintf('Spatial info: %.3f bits/spike, p=%.4f (nShuffles=%d), detected %d valid fields.\n', info_obs, pval, nS, numel(valid_fields));
end

end

% Small helper (keeps single-threaded for now)
function forloop(n, f)
for i=1:n, f(i); end
end
