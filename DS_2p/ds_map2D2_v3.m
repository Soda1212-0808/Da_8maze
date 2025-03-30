clear all
% 定义包含CSV文件的文件夹路径
Path = 'G:\CA3_rawdata\CA3_2p\data';    % 设置数据存放的文件夹路径
% animals={'1464'};
animals={'1646','1306','1307','1309','1311','1312','1974','1976'};

for curr_animal=1:length(animals)
    animal=animals{curr_animal};

    newfolderName = 'buffer_image';
    if exist(fullfile(Path,animal,newfolderName), 'dir') ~= 7
        mkdir(fullfile(Path,animal,newfolderName));
        disp(['Folder "', newfolderName, '" created.']);
    end
    bufferfolderName = 'bufferFile';

    contents = dir(fullfile(Path ,animal));
    % 获取所有子文件夹的名称
    recording_files = {contents([contents.isdir] & ~ismember({contents.name}, {'.', '..'}) & ...
        ~cellfun(@isempty, regexp({contents.name}, '^\d{4}-\d{2}-\d{2}$'))).name}';

    if exist(fullfile(Path, animal ,bufferfolderName), 'dir') ~= 7
        mkdir(fullfile(Path, animal , bufferfolderName));
        disp(['Folder "', bufferfolderName, '" created.']);
    end

    all_data_path=cell(5,1);
    all_data_match=cell(5,1);
    sort_data_path=cell(5,1);
    for curr_day=1:5
        csvFiles = dir(fullfile(Path,animal,recording_files{curr_day},'video_track', '*filtered.csv'))';
        csvFilePaths = fullfile({csvFiles.folder}, {csvFiles.name});
        all_data_path{curr_day} =readtable(csvFilePaths{1});
        all_data_match{curr_day}=load(fullfile(Path,animal,recording_files{curr_day},[recording_files{curr_day},'_video_cell_match.mat']));

        % 根据match 文件 找到cell timepoint对应的 小鼠轨迹
        match_id=cell2mat(all_data_match{curr_day}.animal_match(:,3));
        sort_data_path{curr_day}=all_data_path{curr_day}(match_id(match_id~=0),:);
    end


    %% 保存轨迹 BW1
    bufferfolderName = 'bufferFile';
    X_position=cellfun(@(x)table2array(x(:,2))+100,all_data_path,'UniformOutput',false);
    Y_position=cellfun(@(x)table2array(x(:,3))+100,all_data_path,'UniformOutput',false);

%     ds.draw_trajectory_roi
        load(fullfile(Path, animal ,bufferfolderName,'grab_picture.mat'))


    %% 提取轨迹
    all_data_X_filter_speed=cell(length(X_position),1);
    all_data_Y_filter_speed=cell(length(X_position),1);
    occupancy_time_all=cell(length(X_position),1);
    x_edges=cell(length(X_position),1);
    y_edges=cell(length(X_position),1);
    all_data_speed=cell(length(X_position),1);
    figure('Name',animal,'Position',[50 50 400 800]);
    %设置速度阈值
    speed_threshold=0;
    for curr_day=1:length(X_position)
        X_interp=X_position{curr_day};
        Y_interp=Y_position{curr_day};
        nexttile
        plot(X_interp,Y_interp)
        title(['raw day ' num2str(curr_day) ])

        % 找到 mask 中对应位置为 0 的索引
        mask=BW1{1};
        invalid_indices = mask(sub2ind(size(mask), fix(Y_interp+1), fix(X_interp+1))) == 0;
        % 将对应的 X_interp 和 Y_interp 位置设为 NaN
        X_interp(invalid_indices) = NaN;
        Y_interp(invalid_indices) = NaN;
        nexttile
                plot(X_interp,Y_interp)

        window_size = 5;
        distances_X = abs(X_interp - movmean(X_interp, [window_size window_size], 'omitnan'));
        distances_Y = abs(Y_interp - movmean(Y_interp, [window_size window_size], 'omitnan'));
        % 计算欧几里得距离
        distances = sqrt(distances_X.^2 + distances_Y.^2);
        nexttile;
        plot(distances)
        threshold=100;
        errors = distances > threshold;
        X_interp(errors)=NaN;
        Y_interp(errors)=NaN;

        nexttile
        plot(X_interp,Y_interp);xlim([0 1200]);ylim([0 1200]);
        title(['processed day' num2str(curr_day)])
        speed = [0 ;sqrt(diff(X_interp).^2 + diff(X_interp).^2)];
        nexttile
        plot(speed);
        title('speed')

        X_filter_speed=X_interp; 
%         X_filter_speed(speed<speed_threshold)=NaN;

        Y_filter_speed=Y_interp; 
%         Y_filter_speed(speed<speed_threshold)=NaN;
        
all_data_X_filter_speed{curr_day}=X_filter_speed;
        all_data_Y_filter_speed{curr_day}=Y_filter_speed;
        all_data_speed{curr_day}=speed;

        %%网格分辨率
        bin_size=10;
        % data_path.time=(data_path.scorer)/framerate;
        % position_time= data_path.time;
        frame_sampling_rate=30;
        % 定义网格的边界和分辨率
        % x_edges{curr_day} = min(X_filter_speed):bin_size:max(X_filter_speed);
        % y_edges{curr_day} = min(Y_filter_speed):bin_size:max(Y_filter_speed);
        x_edges{curr_day} = 0:bin_size:1200;
        y_edges{curr_day} = 0:bin_size:1200;

        % 计算占用直方图
        occupancy_map = histcounts2(X_filter_speed, Y_filter_speed, x_edges{curr_day}, y_edges{curr_day});
        % 计算占用时间（假设位置数据的采样率为 position_sampling_rate）
        occupancy_time = occupancy_map * (1 / frame_sampling_rate);

        
        occupancy_time_all{curr_day}=occupancy_time;
    end

    sum(isnan(all_data_X_filter_speed{2}))/length(all_data_X_filter_speed{2})

    %% load Ca2+数据
    for curr_day=1:5
        data_imaging=load(fullfile(Path,animal,recording_files{curr_day},'image_2p',[recording_files{curr_day} '_spikes']));
        all_data_imgaing{curr_day,1}=data_imaging.spikes;
    end
    %% 计算perievent
    %load behavior
    all_data_event=arrayfun(@(curr_day) table2array(readtable(fullfile(Path,animal,recording_files{curr_day},'behavior',[animal '-' recording_files{curr_day} '_repaired778.xlsx']))),1:5,'UniformOutput',false)';

    % 把event的使用video track的帧数 align到ca2+的帧数
    for curr_i=1:length(all_data_match)

        A=double(cell2mat(all_data_match{curr_i}.animal_match(:,3)))';
        B=str2double(all_data_match{curr_i}.animal_match(:,1))';

        idx = find(diff(A) < 0) + 1; % 找到下降点的索引

        if isempty(idx)
            all_data_match{curr_i}.animal_match(:,8) = num2cell(A');
            all_data_match{curr_i}.animal_match(:,7) = num2cell(B');

        else
            % 计算累积调整量，排除初始的0
            adjustments_A = cumsum([0, A(idx-1)]);

            % 利用广播机制创建调整矩阵，并计算总调整量
            delta_A = max((1:numel(A) >= idx(:)) .* adjustments_A(2:end)', [],1);
            % 应用调整
            all_data_match{curr_i}.animal_match(:,8) = num2cell((A + delta_A)');

            adjustments_B = cumsum([0, B(idx-1)]);
            % 利用广播机制创建调整矩阵，并计算总调整量
            delta_B = max((1:numel(B) >= idx(:)) .* adjustments_B(2:end)',[], 1);
            % 应用调整
            all_data_match{curr_i}.animal_match(:,7) = num2cell((B + delta_B)');
        end
    end

    [~,idxx]=cellfun(@(A)   unique(double(cell2mat(A.animal_match(:, 8))), 'first')  ,all_data_match,'UniformOutput',false);
    all_data_frame = cellfun(@(A,B,C) [B(:,1:3) round(interp1(double(cell2mat(A.animal_match(C, 8))), cell2mat(A.animal_match(C, 7)),B(:,4:11), 'linear', 'extrap')) ],...
        all_data_match, all_data_event,idxx,'UniformOutput',false);


    bbX=cellfun(@(x,y)   x(cell2mat(y.animal_match(:,8))),all_data_X_filter_speed,all_data_match,'UniformOutput',false);
    bbY=cellfun(@(x,y)   x(cell2mat(y.animal_match(:,8))),all_data_Y_filter_speed,all_data_match,'UniformOutput',false);

        % 绘制每一次trial sample& choice run的轨迹
        for i=1:5
            figure
            set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
%             bbX=all_data_X_filter_speed{i}(cell2mat(all_data_match{i}.animal_match(:,8)));
%             bbY=all_data_Y_filter_speed{i}(cell2mat(all_data_match{i}.animal_match(:,8)));

            for curr_trial=1:size(all_data_frame{i},1)
                x_b{i}{curr_trial}=bbX{i}(all_data_frame{i}(curr_trial,8):all_data_frame{i}(curr_trial,10));
                y_b{i}{curr_trial}=bbY{i}(all_data_frame{i}(curr_trial,8):all_data_frame{i}(curr_trial,10));
    
                %  x_b{i}{curr_trial}=X_position{i}(all_data_frame{i}(curr_trial,4):all_data_frame{i}(curr_trial,6));
                % y_b{i}{curr_trial}=Y_position{i}(all_data_frame{i}(curr_trial,4):all_data_frame{i}(curr_trial,6));
    
    
                nexttile
                imagesc(BW1{1, 1}  )
                hold on
                plot(x_b{i}{curr_trial},y_b{i}{curr_trial})
                title(num2str(curr_trial))
            end
            sgtitle([animal '-' num2str(i-1) 'sample begin2sample reward'])
            saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal '-' num2str(i-1) 'sample begin2sample reward.jpg']),'jpg')
    
        end

        for curr_day=1:5
            spike_time=all_data_imgaing{curr_day};
            event_time=all_data_frame{curr_day}(:,4:11);
            event_type=all_data_frame{curr_day}(:,2);
            [psth_average(:,curr_day),psth_single_trial(:,curr_day),psth_single_idx(:,curr_day),used_t]=...
                ds.psth(spike_time,event_time,event_type,'smoothing',10,'norm_window',[0 2]);


            SD(curr_day,:)=std(smoothdata(spike_time,2,'gaussian',10),[],2);
        end

    %%


        for curr_cell=1: size(all_data_imgaing{curr_day},1)

            figure('Position',[50 10 1800 1000]);
            tt = tiledlayout(1,length(all_data_frame)*2,'TileSpacing','tight');
            sgtitle(num2str(curr_cell))
            for curr_day=1:length(all_data_frame)
                t_animal = tiledlayout(tt,9,1);
                t_animal.Layout.Tile = curr_day;

                ax1=nexttile(t_animal);

                % 绘制热区图
                spike_times=all_data_imgaing{curr_day}(curr_cell,:)';
                % 获取每个发放事件对应的位置索引
                sample_phase_segments = arrayfun(@(i) all_data_frame{curr_day}(i,4):all_data_frame{curr_day}(i,6), 1:size(all_data_frame{curr_day},1), 'UniformOutput', false);

                spike_x{curr_day} = bbX{curr_day}([sample_phase_segments{:}]');
                spike_y{curr_day} = bbY{curr_day}([sample_phase_segments{:}]');
                spike_rate{curr_day}=spike_times([sample_phase_segments{:}]');
                % 计算发放直方图

                x_idx = discretize(spike_x{curr_day}, x_edges{curr_day});
                y_idx = discretize(spike_y{curr_day}, y_edges{curr_day});

                valid_idx = ~(isnan(x_idx) | isnan(y_idx)); % 仅保留有效索引
                x_idx = x_idx(valid_idx);
                y_idx = y_idx(valid_idx);
                A_valid = spike_rate{curr_day}(valid_idx);

                spike_map{curr_day} = accumarray(...
                    [x_idx, y_idx], ... % 获取每个时间点的格子索引
                    A_valid, ... % 使用放电次数作为权重
                    [length(x_edges{curr_day})-1, length(y_edges{curr_day})-1], ... % 维度
                    @sum, single(0)); % 默认值为 0


%               spike_map{curr_day} = histcounts2(spike_x{curr_day}, spike_y{curr_day}, x_edges{curr_day}, y_edges{curr_day});
                % 计算发放速率地图
                rate_map{curr_day} = spike_map{curr_day} ./ occupancy_time_all{curr_day};
                rate_map{curr_day}(isnan(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0
                rate_map{curr_day}(isinf(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0
                smooth_sigma = 2; % 根据需要调整
                % 对发放速率地图进行高斯平滑
                rate_map_smoothed{curr_day} = imgaussfilt(rate_map{curr_day}, smooth_sigma);
                % 计算平均发放速率
                mean_rate{curr_day} = mean(rate_map{curr_day}(rate_map{curr_day}>0));
                % 识别发放场
                threshold = 1.5 * mean_rate{curr_day}; % 设置阈值为平均发放速率的两倍
                firing_field{curr_day} = rate_map{curr_day} > threshold;
                %%是否反转图像
                %     smoothed_rate_map=flipud(smoothed_rate_map);
                rate_map_smoothed{curr_day}=rate_map_smoothed{curr_day};

                imagesc(x_edges{curr_day}, y_edges{curr_day}, rate_map_smoothed{curr_day});
                axis image off;
                if max(rate_map_smoothed{1})>0;
                    clim([0 max(rate_map_smoothed{1}(:))])
                end
                xlim([0 1200]);     ylim([0 1200]);
                colormap(ax1,'jet')


                for curr_event=1:4

                    aa=nexttile(t_animal)


                %%计算 放电p
                    [idx_type,idx_seq]=sort(psth_single_idx{curr_event,curr_day});
                    heatmap1= sortrows([psth_single_idx{curr_event,curr_day},...
                        permute(psth_single_trial{curr_event,curr_day}(curr_cell,:,:), [2 3 1])],1);

                    heatmap1_baseline=sortrows([psth_single_idx{8,curr_day},...
                        permute(psth_single_trial{8,curr_day}(curr_cell,:,51:75), [2 3 1])],1);
                    heatmap1_baseline1=sum(heatmap1_baseline(:,2:end),2);
                    n_shuff = 1000;  % 置换检验次数

                    shifts=[51 :75];
                    curr_shift=1
                    heatmap_shift=sum(heatmap1(:,shifts),2)
%                     figure;
%                     plot(heatmap1_baseline1);
%                     hold on
%                     plot(heatmap_shift)

                    algine_heatmap=[heatmap1_baseline1 heatmap_shift];

                    % === 进行置换检验 ===
                     types_all=unique(heatmap1(:,1), 'first');
                    event_response_p=zeros(length(types_all),1);
                    for curr_arm =1:length(types_all)
                        event_response = squeeze(mean(diff(algine_heatmap(idx_type==types_all(curr_arm),:),[],2),1));

                        event_response_shuff11 = cell2mat(arrayfun(@(shuff) ...
                            squeeze(mean(diff(ap.shake(algine_heatmap(idx_type==types_all(curr_arm),:),2),[],2),1)), ...
                            1:n_shuff, 'uni', false));

                        event_response_rank = tiedrank(horzcat(event_response, event_response_shuff11)')';
                        event_response_p(curr_arm)= event_response_rank(:,1) ./ (n_shuff + 1)>0.95;
                    end


                    imagesc(used_t,[],heatmap1(:,2:end))

                    hold on

                    mark_dot=find((heatmap1(:,1) .* event_response_p(heatmap1(:,1)))~=0);
                    scatter(zeros(length(mark_dot),1),mark_dot,10,'red','filled')


                    [~, first_indices] = unique(heatmap1(:,1), 'first');
                    first_indices(event_response_p==1)

                    yline([first_indices-0.5])
                    clim([ 0 ...
                        20*max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
                   
                    colormap(aa,ap.colormap('WR'))
                    ax=nexttile(t_animal);

                    if curr_day==1 ||((curr_day==4||curr_day==5) &(curr_event==1||curr_event==2||curr_event==3||curr_event==4))
                        custom_colors = [0 0 0; 0.5 0.5 0.5];
                    elseif curr_day==2||curr_day==3;
                        custom_colors = [0 0 1; 0 0 0; 1 0 0;1 0.5 0.5; 0.5 0.5 0.5;0.5 0.5 1];
                    elseif (curr_day==4||curr_day==5) &(curr_event==7||curr_event==8||curr_event==5||curr_event==6)
                        custom_colors = [0 0 0; 1 0 0; 0 1 0; 1 0.5 0.5; 0.5 0.5 0.5];
                    end
                    ax.ColorOrder = custom_colors; % 设置颜色顺序
                    ax.NextPlot = 'replacechildren'; % 确保颜色按顺序应用

                    h= plot(used_t,psth_average{curr_event,curr_day}(:,:,curr_cell));
                    hold on
                    yline(SD(curr_day,curr_cell))
                    xlim([ -2 2])
ylim([  min(cellfun(@(x) min(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all') ...
                        max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
                                       %                 title([all_event(curr_event).name num2str(curr_day)])
                
                end

            end
            
            for curr_day=1:length(all_data_frame)
               t_animal = tiledlayout(tt,9,1);
                t_animal.Layout.Tile = 5+curr_day;

                ax1=nexttile(t_animal);

                % 绘制热区图
                spike_times=all_data_imgaing{curr_day}(curr_cell,:)';
                % 获取每个发放事件对应的位置索引
                sample_phase_segments = arrayfun(@(i) all_data_frame{curr_day}(i,8):all_data_frame{curr_day}(i,10), 1:size(all_data_frame{curr_day},1), 'UniformOutput', false);

              
                spike_x{curr_day} = bbX{curr_day}([sample_phase_segments{:}]');
                spike_y{curr_day} = bbY{curr_day}([sample_phase_segments{:}]');
                spike_rate{curr_day}=spike_times([sample_phase_segments{:}]');
                % 计算发放直方图

                x_idx = discretize(spike_x{curr_day}, x_edges{curr_day});
                y_idx = discretize(spike_y{curr_day}, y_edges{curr_day});

                valid_idx = ~(isnan(x_idx) | isnan(y_idx)); % 仅保留有效索引
                x_idx = x_idx(valid_idx);
                y_idx = y_idx(valid_idx);
                A_valid = spike_rate{curr_day}(valid_idx);

                spike_map{curr_day} = accumarray(...
                    [x_idx, y_idx], ... % 获取每个时间点的格子索引
                    A_valid, ... % 使用放电次数作为权重
                    [length(x_edges{curr_day})-1, length(y_edges{curr_day})-1], ... % 维度
                    @sum, single(0)); % 默认值为 0



                % 计算发放直方图
%                 spike_map{curr_day} = histcounts2(spike_x{curr_day}, spike_y{curr_day}, x_edges{curr_day}, y_edges{curr_day});
                % 计算发放速率地图
                rate_map{curr_day} = spike_map{curr_day} ./ occupancy_time_all{curr_day};
                rate_map{curr_day}(isnan(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0
                rate_map{curr_day}(isinf(rate_map{curr_day})) = 0; % 将NaN值（由于0占用时间导致的）设为0
                smooth_sigma = 2; % 根据需要调整
                % 对发放速率地图进行高斯平滑
                rate_map_smoothed{curr_day} = imgaussfilt(rate_map{curr_day}, smooth_sigma);
                % 计算平均发放速率
                mean_rate{curr_day} = mean(rate_map{curr_day}(rate_map{curr_day}>0));
                % 识别发放场
                threshold = 1.5 * mean_rate{curr_day}; % 设置阈值为平均发放速率的两倍
                firing_field{curr_day} = rate_map{curr_day} > threshold;
                %%是否反转图像
                %     smoothed_rate_map=flipud(smoothed_rate_map);
                rate_map_smoothed{curr_day}=rate_map_smoothed{curr_day};

                imagesc(x_edges{curr_day}, y_edges{curr_day}, rate_map_smoothed{curr_day});
                axis image off;
                if max(rate_map_smoothed{1})>0;
                    clim([0 max(rate_map_smoothed{1}(:))])
                end
                xlim([0 1200]);     ylim([0 1200]);
                colormap(ax1,'jet')


                for curr_event=5:8

                    
                    
                    aa=nexttile(t_animal)

 %%计算 放电p
                    [idx_type,idx_seq]=sort(psth_single_idx{curr_event,curr_day});
                    heatmap1= sortrows([psth_single_idx{curr_event,curr_day},...
                        permute(psth_single_trial{curr_event,curr_day}(curr_cell,:,:), [2 3 1])],1);

                    heatmap1_baseline=sortrows([psth_single_idx{8,curr_day},...
                        permute(psth_single_trial{8,curr_day}(curr_cell,:,51:75), [2 3 1])],1);
                    heatmap1_baseline1=sum(heatmap1_baseline(:,2:end),2);
                    n_shuff = 1000;  % 置换检验次数

                    shifts=[51 :75];
                    curr_shift=1
                    heatmap_shift=sum(heatmap1(:,shifts),2)
%                     figure;
%                     plot(heatmap1_baseline1);
%                     hold on
%                     plot(heatmap_shift)

                    algine_heatmap=[heatmap1_baseline1 heatmap_shift];

                    % === 进行置换检验 ===
                     types_all=unique(heatmap1(:,1), 'first');
                    event_response_p=zeros(length(types_all),1);
                    for curr_arm =1:length(types_all)
                        event_response = squeeze(mean(diff(algine_heatmap(idx_type==types_all(curr_arm),:),[],2),1));

                        event_response_shuff11 = cell2mat(arrayfun(@(shuff) ...
                            squeeze(mean(diff(ap.shake(algine_heatmap(idx_type==types_all(curr_arm),:),2),[],2),1)), ...
                            1:n_shuff, 'uni', false));

                        event_response_rank = tiedrank(horzcat(event_response, event_response_shuff11)')';
                        event_response_p(curr_arm)= event_response_rank(:,1) ./ (n_shuff + 1)>0.95;
                    end



                    imagesc(used_t,[],heatmap1(:,2:end))

                        hold on

                    mark_dot=find((heatmap1(:,1) .* event_response_p(heatmap1(:,1)))~=0);
                    scatter(zeros(length(mark_dot),1),mark_dot,10,'red','filled')




                    [~, first_indices] = unique(heatmap1(:,1), 'first');
                    hold on
                    yline([first_indices-0.5])
                    clim([  0 ...
                        20*max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
                   
                    colormap(aa,ap.colormap('WR'))
                    ax=nexttile(t_animal);

                    if curr_day==1 ||((curr_day==4||curr_day==5) &(curr_event==1||curr_event==2||curr_event==3||curr_event==4))
                        custom_colors = [0 0 0; 0.5 0.5 0.5];
                    elseif curr_day==2||curr_day==3;
                        custom_colors = [0 0 1; 0 0 0; 1 0 0;1 0.5 0.5; 0.5 0.5 0.5;0.5 0.5 1];
                    elseif (curr_day==4||curr_day==5) &(curr_event==7||curr_event==8||curr_event==5||curr_event==6)
                        custom_colors = [0 0 0; 1 0 0; 0 1 0; 1 0.5 0.5; 0.5 0.5 0.5];
                    end
                    ax.ColorOrder = custom_colors; % 设置颜色顺序
                    ax.NextPlot = 'replacechildren'; % 确保颜色按顺序应用

                    h= plot(used_t,psth_average{curr_event,curr_day}(:,:,curr_cell));
                    hold on
                    yline(SD(curr_day,curr_cell))
                    ylim([  min(cellfun(@(x) min(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all') ...
                        max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
                    xlim([-2 2])
                    %                 title([all_event(curr_event).name num2str(curr_day)])
                
                end

            end
           
            drawnow
           saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal 'neuron' num2str(curr_cell) 'perievent.jpg']),'jpg')

        end


    close all








        %% 计算responsive neurons
% === 参数设定 ===
n_shuff = 1000;  % 置换检验次数
response_t_list = {[-5 0], [0 5], [5 10]};  % 多个 response 时间窗口
baseline_t=[5 10];
n_windows = numel(response_t_list);
alpha_corrected = 0.05 / n_windows;  % Bonferroni 校正阈值

% === 初始化存储变量 ===
C = cell(1,5); % 第一层 (1×5 cell)
for i = 1:5
    C{i} = cell(1,8); % 第二层 (1×8 cell)
    for j = 1:8
        C{i}{j} = cell(1,3); % 第三层 (1×3 cell)
        for k = 1:3
            C{i}{j}{k} = cell(1,7); % 第四层 (1×7 cell)
            for l = 1:7
                C{i}{j}{k}{l} = nan(144,1); % 第五层 (144×1 NaN 矩阵)
            end
        end
    end
end

event_response_all = C;  
event_response_p_all = C;


curr_day=3;
curr_timepoint=2;
curr_window=2;
curr_arm=3;


% 遍历每一天
for curr_day=1:5
            event_time_all=all_data_frame{curr_day}(:,4:11);
            event_type_all=all_data_frame{curr_day}(:,2);
            event_success_all=all_data_frame{curr_day}(:,3);
            spike_time=all_data_imgaing{curr_day};

% 遍历每个event timepoint
for curr_timepoint=1:8

     if curr_timepoint<=4
    [event_idx,any_event_id_arm]=findgroups(floor(event_type_all(2:end-1) / 10) ) ;
    elseif curr_timepoint>4
    [event_idx,any_event_id_arm]=findgroups(mod(event_type_all(2:end-1) , 10)) ;
     end
     [success_idx,success_type]=findgroups(event_success_all(2:end-1));
     
%      [event_idx,event_type]=findgroups(event_type_all(2:end-1));

% === 遍历不同 response 时间窗口 ===
for curr_window = 1:n_windows
       

for curr_arm=1:7

  if ismember(curr_arm,any_event_id_arm)

    event_used_idx= (event_idx== (find(any_event_id_arm==curr_arm))&success_idx==2);


%   event_used_idx(2:end-1);
    response_t = response_t_list{curr_window};  % 取当前 response 窗口
    response_idx=event_time_all(event_used_idx,curr_timepoint);

    response_bins = response_idx(2:end)+ response_t;

    baseline_idx=event_time_all(event_used_idx,8);
    baseline_bins= baseline_idx(1:end-1)+baseline_t;


    event_bins = [baseline_bins, response_bins];


    result =cat(3,arrayfun(@(unit) cell2mat(arrayfun(@(row) arrayfun(@(col) sum(spike_time(unit,event_bins(row,col):event_bins(row,col+1))), ...
        [1 3]), (1:size(event_bins,1))', 'UniformOutput', false)),1:size(spike_time,1), 'UniformOutput', false));

    event_spikes=cat(3,result{:});


    % 计算真实的 (响应 - 基线) 放电变化
    event_response = squeeze(mean(diff(event_spikes, [], 2), 1));
    event_response_all{curr_day}{curr_timepoint}{curr_window}{curr_arm} = event_response;
    
    % === 进行置换检验 ===
    event_response_shuff = cell2mat(arrayfun(@(shuff) ...
        squeeze(mean(diff(ap.shake(event_spikes,2),[],2),1)), ...
        1:n_shuff, 'uni', false));
    
    event_response_rank = tiedrank(horzcat(event_response, event_response_shuff)')';
    event_response_p_all{curr_day}{curr_timepoint}{curr_window}{curr_arm} = event_response_rank(:,1) ./ (n_shuff + 1);

end
end
end
end
end


event_response_p_all_1=cellfun(@(x) cellfun(@(y) cellfun(@(a) cat(2,a{:}),y,'UniformOutput',false),x,'UniformOutput',false),event_response_p_all,'UniformOutput',false);
event_response_p_all_2=cellfun(@(x) cellfun(@(y) cat(3,y{:}),x,'UniformOutput',false),event_response_p_all_1,'UniformOutput',false);
event_response_p_all_3=cellfun(@(x) cat(4,x{:}),event_response_p_all_2,'UniformOutput',false);
event_response_p_all_4=cat(5,event_response_p_all_3{:});

s=permute(max(event_response_p_all_4,[],3),[1 2 4 5 3]);

for i=1:size(event_response_p_all_5,1)
    
     figure('Position',[50 10 1400 200]);
            tt = tiledlayout(1,5,'TileSpacing','tight');
    for j=1:5
        nexttile
    imagesc(permute(event_response_p_all_5(i,:,:,j)>0.95,[2 3 1 4])')
        title(num2str(j))
        ylabel('timepoint')
        xlabel('arm')

    end
    sgtitle(num2str(i))
    drawnow
    saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal 'neuron' num2str(i) 'response.jpg']),'jpg')

end

% === 选择最显著的时间窗口 ===
[~, best_window_idx] = min(cellfun(@(p) min(p), event_response_p_all));
best_response_t = response_t_list{best_window_idx};

% === 找到所有显著时间窗口（Bonferroni 校正）===
significant_windows = find(cellfun(@(p) any(p < alpha_corrected), event_response_p_all));


    
end

