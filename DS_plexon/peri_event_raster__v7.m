clear all

Path='E:\SDdata\WT\each mice\';
% Path = 'H:\CA3_reprocess\each mice\';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};
newfolderName = 'PSTH';
if exist(fullfile(Path,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end

% set parameters
smooth_window=100;
raster_window = [-2,2];
psth_bin_size = 0.005;
t_bins = raster_window(1):psth_bin_size:raster_window(2);


all_response_idx=cell(length(animals),1);
all_psth=cell(length(animals),1);
all_raster=cell(length(animals),1);
for curr_animal=1:length(animals)
    % 获取文件夹中的所有内容
    animal=animals{curr_animal};
    contents = dir([Path animal]);
    % 获取所有子文件夹的名称
    recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

    for curr_day=1:length(recording_files)
        %处理event数据
        event_file=dir(fullfile(Path, animal ,recording_files{curr_day} ,'*data_m1.csv'));
        data_event=csvread(fullfile(event_file.folder , event_file.name));

        events_name= unique (data_event(:,1));
        event_times = arrayfun(@(event) arrayfun(@(id) sort(data_event(data_event(:,1) == id, event)),...
            events_name, 'UniformOutput', false),...
            3:10, 'UniformOutput', false);

        neuron_files=dir(fullfile(Path ,animal , recording_files{curr_day} ,'*.t'));
        spikes_all = arrayfun(@(f) readmclusttfile(fullfile(f.folder, f.name))'/10000, ...
            neuron_files, 'UniformOutput', false);

        templates_all = cell2mat(arrayfun(@(idx, n) repmat(idx, n, 1),...
            (1:numel(cellfun(@numel, spikes_all)))', cellfun(@numel, spikes_all), 'UniformOutput', false));
        [spike_times_timelite, sort_idx] = sort(cell2mat(spikes_all));
        spike_templates = templates_all(sort_idx);


       
        [psth_smooth,raster,~]=...
            cellfun(@(y) cellfun(@(x) ap.psth(spike_times_timelite,x,spike_templates,'smoothing',smooth_window,...
            'window',raster_window,'bin_size',psth_bin_size),y,'UniformOutput',false)...
            ,event_times,'UniformOutput',false);



        psth_smooth_norm=cellfun(@(y) cellfun(@(x) normalize(x,2),y,'UniformOutput',false),...
            psth_smooth,'UniformOutput',false);


        response_unit= cellfun(@(y) cellfun(@(x) x>2,y,'UniformOutput',false),...
            psth_smooth_norm,'UniformOutput',false)

        % responsive neurons
        min_len = 20;

        % 卷积检测每行中是否存在长度为4的连续1
        is_match = cellfun(@(y)  cellfun(@(x) conv2(x(:,find(t_bins>-0.5&t_bins<1)), ones(1, min_len), 'valid') == min_len,...
            y,'UniformOutput',false),...
            response_unit,'UniformOutput',false);  % 结果为逻辑矩阵

        % 对每一行找第一个为 true 的位置（用 max 找 true 的列索引）
        [pos_val, response_idx] = cellfun(@(y)  cellfun(@(x) max(x, [], 2),y,'UniformOutput',false),...
            is_match,'UniformOutput',false);  % pos_val 是是否找到，pos_idx 是起始列
    
        figure
        tiledlayout(6, 8);

        for curr_event=1:6
            for curr_timepoint=1:8

                nexttile
                temp_data=psth_smooth_norm{curr_timepoint}{curr_event}(response_idx{curr_timepoint}{curr_event}>1,:);


                [a,b]=sort(response_idx{curr_timepoint}{curr_event}(response_idx{curr_timepoint}{curr_event}>1));

                imagesc(t_bins ,[],temp_data(b,:))

                clim([2 2.2])
                colormap(ap.colormap('WK'))



            end
        end
        sgtitle([animal ' day ' num2str(curr_day)])
        drawnow
        all_response_idx{curr_animal}{curr_day}= ...
            cell2mat(cellfun(@(x) reshape(x, [1, 1, size(x)]), cat(2,response_idx{:}), 'UniformOutput', false));
        all_psth{curr_animal}{curr_day}= ...
            cell2mat(cellfun(@(x) reshape(x, [1, 1, size(x)]), cat(2,psth_smooth_norm{:}), 'UniformOutput', false));
        all_raster{curr_animal}{curr_day}=raster;
    end
end

    



all_psth_curr_day1=cellfun(@(x)  x{1}  ,all_psth,'UniformOutput',false );
all_psth_curr_day1_1=cat(3,all_psth_curr_day1{:});

all_response_idx_curr_day1=cellfun(@(x)  x{1}  ,all_response_idx,'UniformOutput',false );
all_response_idx_curr_day1_1=cat(3,all_response_idx_curr_day1{:});


all_response_idx_curr_day1_align = nan(size(all_response_idx_curr_day1_1));
all_response_idx_curr_day1_align(all_response_idx_curr_day1_1 ~= 1) = t_bins(all_response_idx_curr_day1_1(all_response_idx_curr_day1_1 ~= 1));

%%

figure
tiledlayout(6, 8);

for curr_event=1:6
for curr_timepoint=1:8
  
      nexttile

temp_idx=  permute(all_response_idx_curr_day1_1(curr_event,curr_timepoint,:),[3,2,1]);
if ismember(curr_timepoint,[2 3 6 7])
    selected_idx=temp_idx>100;

else
selected_idx=temp_idx>1;
end


[a,b]=sort(temp_idx(selected_idx));
temp_data=permute(all_psth_curr_day1_1(curr_event,curr_timepoint,...
    selected_idx,:),[3,4,2,1]);



imagesc(t_bins ,[],temp_data(b,:))

      clim([0 2.2])
      colormap(ap.colormap('WK'))






    end
    end






  %% correlation  相关性分析
spike_binning_t = 0.05; % seconds
spike_binning_t_edges = nanmin(spike_times_timelite):spike_binning_t:nanmax(spike_times_timelite);

binned_spikes_depth = zeros(length(spikes_all),length(spike_binning_t_edges)-1);

for curr_cell=1:length(spikes_all)
      binned_spikes_depth(curr_cell,:)=histcounts(spike_times_timelite( ...
        spike_templates==curr_cell),spike_binning_t_edges);
end

mua_corr = corrcoef(binned_spikes_depth');
figure;
imagesc(1:length(spikes_all),1:length(spikes_all),mua_corr);
axis image;
clim([-1,1].*0.5);
colormap(ap.colormap('BWR'))

% 转为距离矩阵
dist = 1 - mua_corr;
Y = squareform(dist);
Z = linkage(Y, 'average');
leafOrder = optimalleaforder(Z, Y);

% 重排
mua_corr_sorted = mua_corr(leafOrder, leafOrder);

% 可视化
figure;
imagesc(mua_corr_sorted);
clim([-1,1].*0.2);
colormap(ap.colormap('BWR'))
 
raster_window = [-2,2];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

  drawnow
