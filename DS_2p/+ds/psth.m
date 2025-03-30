function [psth_average,psth_single_trial,psth_single_idx,t]=psth(spike_time,event_time,event_type,opts)


arguments
spike_time
event_time
event_type
 % Raster window options
    opts.window (2,1) = [-2,2]
    opts.bin_size = 0.04

    % PSTH post-processing options
    opts.smoothing {mustBeNonnegative} = 0
    opts.norm_window (2,1) = [NaN,NaN]
    opts.softnorm {mustBeNonnegative} = 0
end

t = opts.window(1):opts.bin_size:opts.window(2);
frame=t/opts.bin_size; 

any_event= arrayfun(@(x) frame+event_time(:,x) ,1:size(event_time,2),'UniformOutput',false );
any_event1=cat(3,any_event{:});
any_event_id_sample=findgroups(floor(event_type / 10) ) ;
any_event_id_choice=findgroups(mod(event_type , 10)) ;


for curr_timepoint=1:size(event_time,2)

    if  curr_timepoint<=4
        align_id=any_event_id_sample;
    else
        align_id=any_event_id_choice;
    end
    AA = spike_time;

    B = any_event1(2:end-1,:,curr_timepoint);

    C = AA(sub2ind(size(AA), repmat((1:size(AA,1))', 1, size(B,1), size(B,2)), permute(repmat(B, 1, 1, size(AA,1)), [3, 1, 2])));
    aligned_v_avg_sample_1 = arrayfun(@(x) splitapply(@nanmean,permute(C(x,:,:),[2,3,1]),align_id(2:end-1))', 1:size(C,1),'UniformOutput',false);
    aligned_v_avg_{curr_timepoint}=   cat(3,aligned_v_avg_sample_1{:});
    aligned_v_single_trial{curr_timepoint}=permute(C,[3,2,1]);
    aligned_v_single_idx{curr_timepoint}= align_id(2:end-1);

end
psth_average=aligned_v_avg_';
psth_single_trial=aligned_v_single_trial';
psth_single_idx=aligned_v_single_idx';

% Smooth
if opts.smoothing > 0
    psth_average = cellfun(@(x) smoothdata(x,1,'gaussian',opts.smoothing),psth_average,'UniformOutput',false);
%     psth_single_trial=cellfun(@(x) smoothdata(x,1,'gaussian',opts.smoothing),psth_single_trial,'UniformOutput',false);
end

if ~all(isnan(opts.norm_window))
    t_baseline = t >= opts.norm_window(1) & ...
        t <= opts.norm_window(2);
    % (compute baseline as average across all alignments)
    psth_baseline_mean = nanmean(psth_average{8}(t_baseline,:,:),[1,2]);
    psth_baseline_std= std(psth_single_trial{8}(t_baseline,:,:),0,[1,2],'omitnan');
    psth_average =cellfun(@(x) (x - psth_baseline_mean)./psth_baseline_std,psth_average,'UniformOutput',false);
    psth_single_trial=cellfun(@(x) (x - psth_baseline_mean)./psth_baseline_std,psth_single_trial,'UniformOutput',false);
end

