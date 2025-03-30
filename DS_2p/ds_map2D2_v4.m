clear all
% 定义包含CSV文件的文件夹路径
Path = 'G:\CA3_rawdata\CA3_2p\data';    % 设置数据存放的文件夹路径
% animals={'1464'};
animals={'1646','1306','1307','1309','1311','1312','1974','1976'};

for curr_animal=3:length(animals)
                preload_vars = who;

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
%%

%             [B,L] = bwboundaries(BW1{1, 1} );
%            
%   for  curr_cell=1:size(all_data_imgaing{curr_day,1},1)
%     % 绘制每一次trial sample& choice run的轨迹
%     figure('Position',[50 50 1400 1000])
%     
%     for curr_day=1:5
%         for curr_fig=1:7
%             subplot(5, 7, (curr_day-1)*7+curr_fig);
%             hold on;
%             for k = 1:length(B)
%                 boundary = B{k}; % 提取边界坐标
%                 plot(boundary(:,2), boundary(:,1),'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
%             end
%             axis image off;
%         end
% 
%         for curr_trial=1:size(all_data_frame{curr_day},1)
%             x_b{curr_day}{curr_trial}=bbX{curr_day}(all_data_frame{curr_day}(curr_trial,8):all_data_frame{curr_day}(curr_trial,10));
%             y_b{curr_day}{curr_trial}=bbY{curr_day}(all_data_frame{curr_day}(curr_trial,8):all_data_frame{curr_day}(curr_trial,10));
% 
%             spike_times=all_data_imgaing{curr_day}(curr_cell,:)';
%             spike_rate{curr_trial}=spike_times(all_data_frame{curr_day}(curr_trial,8):all_data_frame{curr_day}(curr_trial,10));
%             spike_position_x{curr_trial}=x_b{curr_day}{curr_trial}(spike_rate{curr_trial}>0);
%             spike_position_y{curr_trial}=y_b{curr_day}{curr_trial}(spike_rate{curr_trial}>0);
%             spike_intensity{curr_trial}=spike_rate{curr_trial}(spike_rate{curr_trial}>0);
%             
%             % 设定强度阈值
%             threshold = max(spike_intensity{curr_trial});
%             % 归一化强度（限制最大值为100）
%             normalized_intensity = spike_intensity{curr_trial};
%             normalized_intensity(spike_intensity{curr_trial} > threshold) = threshold; % 超过100的设为100
%             normalized_intensity = normalized_intensity / threshold; % 归一化到 [0,1]
% 
%             
%             figure_idx= mod(all_data_frame{curr_day}(curr_trial,2),10);
%             subplot(5,7,(curr_day-1)*7+figure_idx)
% 
%             scatter(spike_position_x{curr_trial},spike_position_y{curr_trial},10,normalized_intensity,'filled', 'MarkerFaceAlpha', 0.5)
%             colormap( flipud(gray)); % 仅对 scatter 设置灰度 colormap
%             
%             caxis( [0 1]); % 颜色范围 0（白）到 1（黑）
% %             title(num2str(curr_day))
% 
%             
% 
%         end
%         
% %         sgtitle([animal ' cell ' num2str(curr_cell) ' day-' num2str(i-1) 'choice begin2choice reward'])
% %         saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal ' cell' num2str(curr_cell) ' day-' num2str(i-1) 'choice begin2choice reward.jpg']),'jpg')
% 
%     end
%      sgtitle([animal ' cell ' num2str(curr_cell)  'choice begin2choice reward'])
%      drawnow
% 
%         saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal ' cell' num2str(curr_cell)  'choice begin2choice reward.jpg']),'jpg')
% 
%     close all
%   end
%                 

%%
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
            %             mean_rate{curr_day} = mean(rate_map{curr_day}(rate_map{curr_day}>0));
            % 识别发放场
            %             threshold = 1.5 * mean_rate{curr_day}; % 设置阈值为平均发放速率的两倍
            %             firing_field{curr_day} = rate_map{curr_day} > threshold;
            %%是否反转图像
          %   smoothed_rate_map=flipud(smoothed_rate_map);

            imagesc(x_edges{curr_day}, y_edges{curr_day}, rate_map_smoothed{curr_day});
            axis image off;
            if max(rate_map_smoothed{1})>0;
                clim([0 max(rate_map_smoothed{1}(:))])
            end
            xlim([0 1200]);     ylim([0 1200]);
            colormap(ax1,'jet')


            for curr_event=1:4
                aa=nexttile(t_animal);
                heatmap1= sortrows([psth_single_idx{curr_event,curr_day},...
                    psth_single_trial{curr_event,curr_day}(:,:,curr_cell)'],1);
                imagesc(used_t,[],heatmap1(:,2:end))
                hold on

                [~, first_indices] = unique(heatmap1(:,1), 'first');

                yline([first_indices-0.5])
          
                clim([0 5])
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

                % nexttile
                % imagesc(psth_average{curr_event,curr_day}(:,:,curr_cell)>2*SD(curr_day,curr_cell))
                %


                hold on
%                 yline(2*SD(curr_day,curr_cell),'Color','r')
                xlim([ -2 2])
%                 ylim([  min(cellfun(@(x) min(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all') ...
%                     max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
ylim([-0.2 4])
yline(1.96)
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


                heatmap1= sortrows([psth_single_idx{curr_event,curr_day},...
                    psth_single_trial{curr_event,curr_day}(:,:,curr_cell)'],1);


                imagesc(used_t,[],heatmap1(:,2:end))



                [~, first_indices] = unique(heatmap1(:,1), 'first');
                hold on
                yline([first_indices-0.5])
%                 clim([ 0 ...
%                     0.5*max(heatmap1(2:end),[],'all')])
clim([ 0 5])
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
                %                 yline(2*SD(curr_day,curr_cell),'Color','r')
                %                 ylim([  min(cellfun(@(x) min(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all') ...
                %                     max(cellfun(@(x)  max(x(:,:,curr_cell),[],"all"),psth_average,'UniformOutput',true),[],'all')])
                xlim([-2 2])
                ylim([-0.2 4])
                yline(1.96)

            end

        end

        drawnow
        saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal 'neuron' num2str(curr_cell) 'perievent.jpg']),'jpg')

    end


    close all

            clearvars('-except',preload_vars{:});



end

