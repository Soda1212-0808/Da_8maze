function [spike_x, spike_y, spike_trial_idx] = map_spikes_to_positions(...
    spike_times, position_timelite, position_x, position_y, inIntervals, frame_rate, varargin)
% MAP_SPIKES_TO_POSITIONS_BY_INTERVALS
% Map spike times -> positions using position vectors split by intervals (inIntervals is T x n logical).
%
% USAGE
%   [sx, sy, strial] = map_spikes_to_positions_by_intervals(spike_times, position_timelist, ...
%           position_x, position_y, inIntervals, frame_rate)
%
% INPUTS
%   spike_times        - Sx1 vector (seconds)
%   position_timelist  - Tx1 vector (seconds)
%   position_x, position_y - Tx1 vectors
%   inIntervals        - T x n logical matrix. Each column contains a (preferably) contiguous block of true
%                        indicating that column's trial/time-segment.
%   frame_rate         - scalar, frames per second (used only if you need to extend end-time by 1/frame_rate)
%
% Name-Value options:
%   'method'           - interp method: 'linear' (default), 'spline', 'nearest'
%   'allow_extrap'     - logical, default false. If true, interp1 extrapolates; otherwise out-of-range spikes -> NaN
%   'fill_nan'         - logical, default true. If true, fill NaNs inside each trial trajectory by local interp1 before mapping.
%   'verbose'          - logical, default false
%
% OUTPUTS
%   spike_x            - Sx1 vector, x positions for each spike (NaN if not mapped)
%   spike_y            - Sx1 vector, y positions
%   spike_trial_idx    - Sx1 vector, trial index (1..n) the spike was mapped to, 0 if none
%
% NOTES
%   - For each column (trial) we take the indices idx = find(inIntervals(:,trial));
%     times = position_timelist(idx); x = position_x(idx); y = position_y(idx).
%     If 'fill_nan' true we fill internal NaNs in x/y via interp1 on valid samples.
%   - We then find spikes with spike_times in [times(1), times(end)] and interp1 within that trial.
%   - This avoids interpolation across gaps between trials.

% parse name-value
p = inputParser;
addParameter(p,'method','linear', @(s) ismember(s,{'linear','spline','nearest'}));
addParameter(p,'allow_extrap',false, @(v)islogical(v) || (isnumeric(v)&&isscalar(v)));
addParameter(p,'fill_nan', true, @(v)islogical(v) || (isnumeric(v)&&isscalar(v)));
addParameter(p,'verbose', false, @(v)islogical(v) || (isnumeric(v)&&isscalar(v)));
parse(p, varargin{:});
method = p.Results.method;
allow_extrap = logical(p.Results.allow_extrap);
fill_nan = logical(p.Results.fill_nan);
verbose = logical(p.Results.verbose);

% ensure column vectors
spike_times = spike_times(:);
position_timelite = position_timelite(:);
position_x = position_x(:);
position_y = position_y(:);

S = numel(spike_times);
T = numel(position_timelite);
if numel(position_x) ~= T || numel(position_y) ~= T
    error('position_timelist, position_x and position_y must have same length T.');
end
if ~islogical(inIntervals)
    inIntervals = inIntervals ~= 0;
end
if size(inIntervals,1) ~= T
    error('inIntervals must have T rows (same length as position_timelist).');
end

% prepare outputs (default NaN / 0)
spike_x = NaN(S,1);
spike_y = NaN(S,1);
spike_trial_idx = zeros(S,1);

% optionally choose extrap val
if allow_extrap
    extrapVal = 'extrap';
else
    extrapVal = NaN;
end

nTrials = size(inIntervals,2);

% For speed, precompute spikes sorted indices (we'll still map by trial)
[spike_times_sorted, sortIdx] = sort(spike_times);
invSortIdx = zeros(S,1); invSortIdx(sortIdx) = 1:S; % to map back
assigned_mask_sorted = false(S,1); % mark which sorted spikes assigned

% For each trial, map spikes within that trial
for tr = 1:nTrials
    idx = find(inIntervals(:,tr));
    if isempty(idx)
        if verbose, warning('Trial %d: empty interval (no indices). Skipping.', tr); end
        continue;
    end
    % ensure contiguous block: find contiguous segments; pick union (we'll handle multiple segments)
    % But here we will accept columns that might have several segments: process each contiguous segment separately
    col = inIntervals(:,tr);
    d = diff([0; col; 0]);
    starts = find(d==1);
    ends = find(d==-1)-1;
    for segi = 1:numel(starts)
        sidx = starts(segi);
        eidx = ends(segi);
        idx_seg = sidx:eidx;
        t_seg = position_timelite(idx_seg);
        x_seg = position_x(idx_seg);
        y_seg = position_y(idx_seg);
        % optionally fill internal NaNs in x_seg/y_seg by local interp (only inside segment)
        if fill_nan
            valid_x = ~isnan(x_seg);
            if sum(valid_x) >= 2
                x_seg = interp1(find(valid_x), x_seg(valid_x), (1:numel(x_seg))', 'linear', NaN);
            end
            valid_y = ~isnan(y_seg);
            if sum(valid_y) >= 2
                y_seg = interp1(find(valid_y), y_seg(valid_y), (1:numel(y_seg))', 'linear', NaN);
            end
            % Note: if fewer than 2 valid points remain, interp1 leaves NaNs -> spikes in that region won't be mapped
        end

        % find spikes in this segment using sorted spike times (fast)
        t0 = t_seg(1);
        t1 = t_seg(end);
        % find sorted spike indices in [t0, t1]
        left = find(spike_times_sorted >= t0, 1, 'first');
        if isempty(left)
            continue;
        end
        right = find(spike_times_sorted <= t1, 1, 'last');
        if isempty(right) || right < left
            continue;
        end
        spikes_idx_sorted = left:right;
        spikes_times_in = spike_times_sorted(spikes_idx_sorted);

        % interpolate positions for these spikes
        % use interp1 with t_seg as x and x_seg/y_seg as y
        x_vals = interp1(t_seg, x_seg, spikes_times_in, method, extrapVal);
        y_vals = interp1(t_seg, y_seg, spikes_times_in, method, extrapVal);

        % keep only those with non-NaN positions (interp may produce NaN if insufficient data)
        valid_spikes = ~isnan(x_vals) & ~isnan(y_vals);
        if ~any(valid_spikes)
            continue;
        end

        % map back to original spike ordering
        spikes_sorted_valid = spikes_idx_sorted(valid_spikes);         % positions in sorted array
        spikes_orig_idx = sortIdx(spikes_sorted_valid);               % indices into original spike_times
        % assign outputs
        spike_x(spikes_orig_idx) = x_vals(valid_spikes);
        spike_y(spikes_orig_idx) = y_vals(valid_spikes);
        spike_trial_idx(spikes_orig_idx) = tr;
        assigned_mask_sorted(spikes_sorted_valid) = true;
    end
end

% optionally report unassigned spikes
if verbose
    nAssigned = sum(spike_trial_idx>0);
    nUnassigned = S - nAssigned;
    fprintf('Mapped %d / %d spikes to positions; %d unassigned (outside intervals or unmappable).\n', nAssigned, S, nUnassigned);
end

end
