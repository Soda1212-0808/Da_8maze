function [psth_average, psth_single_trial, t] = psth_2p_v2(spike_time, event_time, opts)
% PSTH  compute peri-event time histogram (PETH) / average activity for frame-aligned calcium traces
%
% USAGE:
%   [psth_avg, psth_trials, t] = psth(spike_time, event_time)
%   [psth_avg, psth_trials, t] = psth(spike_time, event_time, 'frame_rate', 15, 'window', [-1 1], ...)
%
% INPUTS:
%   spike_time - T x 1 vector: per-frame activity (deconvolved amplitude or dF/F per frame)
%   event_time - n x 1 vector: event frame indices (1-based). Each entry is an integer frame index.
%
% OPTIONS (name–value via opts):
%   opts.frame_rate   - frames per second (Hz). Default 10.
%   opts.window       - 1x2 in seconds (default [-1 1])
%   opts.smoothing    - smoothing window for smoothdata (in #samples/frames). Default 3
%   opts.do_smooth    - whether to smooth (default true)
%   opts.norm_window  - 1x2 in seconds for baseline normalization (default [-1 0])
%
% OUTPUTS:
%   psth_average      - nBins x 1 average across events (z-scored if norm applied)
%   psth_single_trial - nBins x nEvents matrix, each column is one aligned trial
%   t                 - time vector relative to event (seconds). Step = 1/frame_rate

arguments
    spike_time (:,1) double
    event_time (:,1) double
    opts.frame_rate (1,1) double = 10
    opts.window (1,2) double = [-1 1]
    opts.smoothing (1,1) double {mustBeNonnegative} = 3
    opts.do_smooth (1,1) logical = true
    opts.norm_window (1,2) double = [-1 0]
end

% ---- basic checks ----
T = numel(spike_time);
if ~isvector(spike_time) || size(spike_time,2)~=1
    error('spike_time must be a T-by-1 vector.');
end

% ensure integer frame indices, drop out-of-range events
event_time = round(event_time(:));
valid_mask_events = (event_time >= 1) & (event_time <= T);
event_time = event_time(valid_mask_events);
nEvents = numel(event_time);
if nEvents == 0
    error('No valid events within [1, %d] frames after filtering.', T);
end

frame_rate = opts.frame_rate;

% ---- time axis: step = 1/frame_rate (one entry per frame relative to event) ----
t = opts.window(1) : (1/frame_rate) : opts.window(2);
nBins = numel(t);

% offsets in frames (integers)
frameOffset = round(t * frame_rate);  % e.g., for frame_rate=10, t step 0.1s -> offsets -10:-9:...:10

% allocate outputs
psth_single_trial = NaN(nBins, nEvents);

% ---- build PSTH per trial by sampling spike_time at event_frame + offset ----
for k = 1:nEvents
    eframe = event_time(k);
    idx = eframe + frameOffset;            % indices for each bin (may be out of range)
    validIdx = (idx >= 1) & (idx <= T);
    if any(validIdx)
        psth_single_trial(validIdx, k) = spike_time(idx(validIdx));
    end
    % out-of-range bins remain NaN
end

% ---- average across trials (ignoring NaNs) ----
psth_average = nanmean(psth_single_trial, 2);

% ---- smoothing (along time dimension) ----
if opts.do_smooth && opts.smoothing > 0
    w = max(1, round(opts.smoothing)); % smoothing window in frames
    % smooth average
    psth_average = smoothdata(psth_average, 1, 'gaussian', w);
    % (optionally smooth single trials too — uncomment if desired)
    % for k = 1:nEvents
    %     psth_single_trial(:,k) = smoothdata(psth_single_trial(:,k), 1, 'gaussian', w);
    % end
end

% ---- baseline normalization (z-scoring) across baseline window if requested ----
if ~all(isnan(opts.norm_window))
    baseline_mask = (t >= opts.norm_window(1)) & (t <= opts.norm_window(2));
    if any(baseline_mask)
        baseline_vals = psth_single_trial(baseline_mask, :);
        baseline_vec = baseline_vals(:);
        baseline_vec = baseline_vec(~isnan(baseline_vec));
        if ~isempty(baseline_vec)
            mu = mean(baseline_vec);
            sigma = std(baseline_vec);
            if sigma == 0
                % subtract mean only to avoid division by zero
                psth_average = psth_average - mu;
                psth_single_trial = psth_single_trial - mu;
            else
                psth_average = (psth_average - mu) ./ sigma;
                psth_single_trial = (psth_single_trial - mu) ./ sigma;
            end
        else
            warning('No valid baseline samples across trials; skipping normalization.');
        end
    else
        warning('norm_window does not overlap PSTH time vector; skipping normalization.');
    end
end

end
