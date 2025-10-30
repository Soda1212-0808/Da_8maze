clear all
clc

Path=ds.locations.server_data_path;
load_parts=struct;
load_parts.tetrode=true;
load_parts.video_track=true;
animal='DCA3-20';
% rec_day='2021-06-13';
recordings=ds.find_recordings(animal);
rec_day = recordings(1).day;
ds.load_recoding

time_intervals=cellfun(@(tb,tc) ...
    arrayfun(@(i) position_timelite >= tb(i) & position_timelite <= tc(i), (1:size(tb,1))', ...
    'UniformOutput', false),arm_times{1},arm_times{3},'UniformOutput',false);



%% 计算占用直方图

bin_size=10;
smooth_bin=2;
inIntervals=cell2mat(cat(1, time_intervals{:})');
    colors=num2cell(jet(6), 2);

    % [occupancy_time,x_edges,y_edges]=ds.occupancy(position_re_X,position_re_Y,'timeinterval',inIntervals,'bin_size',10,'fps',30);
    result=cell(length(spikes_all),1);
    for curr_cell=71:length(spikes_all)

        spike_times=spikes_all{curr_cell};

        % [rate_map,  x_edges, y_edges] = ...
        %     ds.compute_rate_map_from_position(position_re_X, position_re_Y, position_timelite, spike_times, inIntervals,bin_size,...
        %     'smooth',true, 'smooth_bin', smooth_bin,'bin_pass',5);

        result{curr_cell}=ds.detect_place_cell(position_re_X, position_re_Y, position_timelite, spike_times, inIntervals,bin_size,...
            'smooth', true, 'smooth_bin', 2, 'bin_pass', 5, 'nShuffles', 500, ...
            'info_thresh', 0.15, 'p_thresh', 0.05, 'min_field_bins', 16, 'peak_rate_min', 1);


        modified_string = strrep(neuron_files(curr_cell).name(1:end-2), '_', '-');
        figure('Name',modified_string,'Position',[50 50 300 800])
        tiledlayout(5,1)
        nexttile
        % imagesc(x_edges, y_edges, result.rate_map_smooth);
        imagesc( result{curr_cell}.rate_map_smooth);

        axis image off;
        clim([0 nanmax(result{curr_cell}.rate_map_smooth(:))])
        colormap("jet")

%         if result{curr_cell}.info>0.1 & result{curr_cell}.p_value<0.05

        hold on
        for curr_field =1:length(result{curr_cell}.fields_valid)
            boundaries = bwboundaries(result{curr_cell}.fields_valid(curr_field).mask);
            plot(boundaries{1}(:,2), boundaries{1}(:,1), 'r', 'LineWidth', 1);
        end

%         end

        formatted_value = sprintf('%.1f', round(spike_freq(curr_cell),1));
        title([modified_string ': ' formatted_value 'Hz'])

        % psth
        % set parameters
        smooth_window=100;
        raster_window = [-2,2];
         psth_bin_size = 0.01;
        % t_bins = raster_window(1):psth_bin_size:raster_window(2);
        % reponse_window=[-0.5 0.5];

        [use_spikes,spike_groups] = ismember(spike_templates,curr_cell);
        [psth_smooth,raster,raster_t_stim]=cellfun(@(y)...
            cellfun(@(x) ap.psth(spike_times_timelite(use_spikes),x,spike_groups(use_spikes),'smoothing',smooth_window,...
            'window',raster_window,'bin_size',psth_bin_size)...
            ,y,'UniformOutput',false),arm_times,'UniformOutput',false);


        for curr_stage=1:2
        nexttile
        hold on
        cellfun(@(x,color) plot(raster_t_stim{1}{1},smoothdata(x,2,'gaussian',10)','linewidth',2,'Color',color),psth_smooth{curr_stage},colors,'UniformOutput',false)
        xlim([-2 2])

        nexttile

        idxMat = repelem((1:numel(raster{curr_stage}))', cellfun(@numel, arm_times{curr_stage}));
        [raster_y,raster_x] = find(cat(1,raster{curr_stage}{:}));
        hold on
        for curr_i=1:length(arm_times{1})

            cur_idx= ismember(raster_y,  find(idxMat==curr_i));

            plot(raster_t_stim{1}{1}(raster_x(cur_idx)),raster_y(cur_idx),'LineStyle','none','Marker','.','Color',colors{curr_i});
        end
        axis off

        end

        drawnow

    end



    result_all=vertcat(result{:});

    find(cell2mat({result_all.info})>0&cell2mat({result_all.p_value})<0.05)

