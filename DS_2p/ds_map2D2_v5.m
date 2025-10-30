clear all
clc

Path=ds.locations.server_data_path;
load_parts=struct;
load_parts.ca_2p=true;
load_parts.video_track=true;
animal='CA3-1307';
% rec_day='2021-06-13';
recordings=ds.find_recordings(animal);
rec_day = recordings(1).day;
ds.load_recoding








%%
position_frame_to_2p=(1:size(animal_match,1))';
position_timelite_to_2p=(1:size(animal_match,1))'/10;
arm_times_to_2p =cellfun(@(y) cellfun(@(x) ...
    round(interp1(double(cell2mat(animal_match(idxx, 8))), cell2mat(animal_match(idxx, 7)),x, 'linear', 'extrap')),...
    y,'UniformOutput',false),arm_times,'UniformOutput',false)
time_intervals=cellfun(@(tb,tc) ...
    arrayfun(@(i) position_frame_to_2p >= tb(i) & position_frame_to_2p <= tc(i), (1:size(tb,1))', ...
    'UniformOutput', false),arm_times_to_2p{1},arm_times_to_2p{3},'UniformOutput',false);
inIntervals=cell2mat(cat(1, time_intervals{:})');


X_to_2p=   position_re_X(cell2mat(animal_match(:,8)));
Y_to_2p=   position_re_Y(cell2mat(animal_match(:,8)));


bin_size=10;
frame_rate =10;

for curr_cell=1:size(ca_2p_spikes,1)
% rate map
 results = ds.detect_place_cell_v2(X_to_2p, Y_to_2p, position_timelite_to_2p, ca_2p_spikes(curr_cell,:), inIntervals, bin_size, ...
     'data_type','ca', 'smooth', true, 'smooth_bin', 2, 'nShuffles', 1000, ...
    'info_thresh', 0.3, 'p_thresh', 0.05, 'peak_rate_min', 0.1,  'verbose', true);

%psth
[psth_average,psth_single_trial,used_t]= ds.psth_2p_v2(ca_2p_spikes(curr_cell,:),arm_times{1, 1}{1, 1}  )


end




