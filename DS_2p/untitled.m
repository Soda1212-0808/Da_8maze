

track_sample_arm=BW1{2};
track_reward_arm= BW1{3}+BW1{4}+BW1{5}+BW1{6}+BW1{7}+BW1{8}+BW1{9};

figure;imagesc(track_reward_arm)

sample_begin =cell(length(all_data_event.all_event_frame),1)
choice_begin =cell(length(all_data_event.all_event_frame),1)
sample_arm =cell(length(all_data_event.all_event_frame),1)
choice_arm =cell(length(all_data_event.all_event_frame),1)

for curr_day=1:length(all_data_event.all_event_frame)

    sample_begin{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);
    choice_begin{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);

    sample_arm{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);
    choice_arm{curr_day}=nan(size(all_data_event.all_event_frame{curr_day},1),1);

    figure;
    % imagesc(track_sample_arm)
    % imagesc(track_reward_arm)

    for curr_trial=1:size(all_data_event.all_event_frame{curr_day},1)

        % sample begin time point
        sample_range=all_data_event.all_event_frame{curr_day}(curr_trial,[1,3]);
        sample_x1=X_filter_speed_all{curr_day}(sample_range(1):sample_range(2));
        sample_y1=Y_filter_speed_all{curr_day}(sample_range(1):sample_range(2));
        sample_in_target_area_1=nan(length(sample_x1),1);

        sample_in_target_area_1(~isnan(sample_x1)) = track_sample_arm(sub2ind(size(track_sample_arm), round(sample_y1(~isnan(sample_x1))), round(sample_x1(~isnan(sample_x1)))));

        if ~isempty(find(sample_in_target_area_1==1,1,'last'))
            sample_begin{curr_day}(curr_trial)=sample_range(1)-1+ find(sample_in_target_area_1==1,1,'last');
        else  sample_begin{curr_day}(curr_trial)=sample_range(1);
        end

        sample_in_target_area_2=nan(length(sample_x1),1);
        sample_in_target_area_2(~isnan(sample_x1)) = track_reward_arm(sub2ind(size(track_reward_arm), round(sample_y1(~isnan(sample_x1))), round(sample_x1(~isnan(sample_x1)))));

        if ~isempty(find(sample_in_target_area_2==1,1,'first'))
            sample_arm{curr_day}(curr_trial)=sample_range(1)-1+ find(sample_in_target_area_2==1,1,'first');
            %     else  sample_arm(curr_trial)=sample_range(1);
        end


        nexttile
        hold on
        imagesc(track_sample_arm+track_reward_arm)
        plot(round(sample_x1), round(sample_y1))
        scatter(sample_x1(find(sample_in_target_area_1==1,1,'last')),sample_y1(find(sample_in_target_area_1==1,1,'last')),'red')
        scatter(sample_x1(find(sample_in_target_area_2==1,1,'first')),sample_y1(find(sample_in_target_area_2==1,1,'first')),'yellow')
        % sample begin time point
        title('sample')


        choice_range=[all_data_event.all_event_frame{curr_day}(curr_trial,4) all_data_event.all_event_frame{curr_day}(curr_trial,5)+100];
        choice_x1=X_filter_speed_all{curr_day}(choice_range(1):choice_range(2));
        choice_y1=Y_filter_speed_all{curr_day}(choice_range(1):choice_range(2));
        choice_in_target_area_1=nan(length(choice_x1),1);
        choice_in_target_area_1(~isnan(choice_x1)) = track_sample_arm(sub2ind(size(track_sample_arm), round(choice_y1(~isnan(choice_x1))), round(choice_x1(~isnan(choice_x1)))));
        if ~isempty(find(choice_in_target_area_1==1,1,'last'))
            choice_begin{curr_day}(curr_trial)=choice_range(1)-1+ find(choice_in_target_area_1==1,1,'last');
        else  choice_begin{curr_day}(curr_trial)=choice_range(1);
        end

        choice_in_target_area_2=nan(length(choice_x1),1);
        choice_in_target_area_2(~isnan(choice_x1)) = track_reward_arm(sub2ind(size(track_reward_arm), round(choice_y1(~isnan(choice_x1))), round(choice_x1(~isnan(choice_x1)))));
        if ~isempty(find(choice_in_target_area_2==1,1,'first'))
            choice_arm{curr_day}(curr_trial)=choice_range(1)-1+ find(choice_in_target_area_2==1,1,'first');
            %     else  choice_arm(curr_trial)=choice_range(1);
        end

        nexttile
        hold on
        imagesc(track_sample_arm+track_reward_arm)
        plot(round(choice_x1), round(choice_y1))
        scatter(choice_x1(find(choice_in_target_area_1==1,1,'last')),choice_y1(find(choice_in_target_area_1==1,1,'last')),'red')
        scatter(choice_x1(find(choice_in_target_area_2==1,1,'first')),choice_y1(find(choice_in_target_area_2==1,1,'first')),'yellow')
        title('choice')
        % sample begin time point

    end


end


for curr_cell=1:size(data_imgaing_all{curr_day},1)
    
    
figure('Position',[50 100 800 400]);
tt = tiledlayout(1,5,'TileSpacing','tight');
sgtitle(num2str(curr_cell))

    for curr_day=1:length(all_data_event.all_event_frame)

        event_range=-50:100;
%sample
        sample_event=event_range+sample_arm{curr_day}(~isnan(sample_arm{curr_day}));
        event_frame=all_data_event.all_event_frame{curr_day}(:,6);
        [idx,order]=sort(event_frame(~isnan(sample_arm{curr_day})));
        [align_id,n]=findgroups(event_frame(~isnan(sample_arm{curr_day})));


        imaging1=data_imgaing_all{curr_day}(curr_cell,:)>50;
        peri_event_sample_imaging{curr_cell}=imaging1(sample_event);

         t_animal = tiledlayout(tt,2,1);
    t_animal.Layout.Tile = curr_day;
%     title(t_animal,num2str(curr_day));
        nexttile(t_animal);

        aligned_v_avg_sample = splitapply(@nanmean,peri_event_sample_imaging{curr_cell},align_id)';

        h= plot(-5:0.1:10, smoothdata(aligned_v_avg_sample,1,'gaussian',20));
        if curr_day==1
            h(1).Color = [0, 0, 0];  % 第一条线，颜色为红色
            h(2).Color = [0.5, 0.5, 0.5];  % 第二条线，颜色为绿色
        elseif  curr_day==2||curr_day==3
            h(1).Color = [0, 0, 1];  % 第一条线，颜色为红色
            h(2).Color = [0, 0, 0];  % 第二条线，颜色为绿色
            h(3).Color = [1, 0, 0];  % 第三条线，颜色为蓝色
            h(4).Color = [1, 0.5, 0.5];  % 第四条线，颜色为黑色
            h(5).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
            h(6).Color = [0.5, 0.5, 1];  % 第五条线，颜色为品红色
        elseif  curr_day==4||curr_day==5
            h(1).Color = [0, 0, 0];  % 第一条线，颜色为红色
            h(2).Color = [0, 0, 0];  % 第二条线，颜色为绿色
            h(3).Color = [0, 0, 0];  % 第三条线，颜色为蓝色
            h(4).Color = [0, 0, 0];  % 第四条线，颜色为黑色
            h(5).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
            h(6).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
            h(7).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
            h(8).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
        end
    
    ylim([0 0.4])
    title(['sample day ' num2str(curr_day)])

%choice
choice_event=event_range+choice_arm{curr_day}(~isnan(choice_arm{curr_day}));
        event_frame=all_data_event.all_event_frame{curr_day}(:,6);
        [idx,order]=sort(event_frame(~isnan(choice_arm{curr_day})));
        [align_id,n]=findgroups(event_frame(~isnan(choice_arm{curr_day})));


        imaging1=data_imgaing_all{curr_day}(curr_cell,:)>50;
        peri_event_choice_imaging{curr_cell}=imaging1(choice_event);
        nexttile(t_animal);
        aligned_v_avg_choice = splitapply(@nanmean,peri_event_choice_imaging{curr_cell},align_id)';
        
        h= plot(-5:0.1:10, smoothdata(aligned_v_avg_choice,1,'gaussian',20));
        if curr_day==1
            h(2).Color = [0, 0, 0];  % 第一条线，颜色为红色
            h(1).Color = [0.5, 0.5, 0.5];  % 第二条线，颜色为绿色
        elseif  curr_day==2||curr_day==3
            h(6).Color = [0, 0, 1];  % 第一条线，颜色为红色
            h(5).Color = [0, 0, 0];  % 第二条线，颜色为绿色
            h(4).Color = [1, 0, 0];  % 第三条线，颜色为蓝色
            h(3).Color = [1, 0.5, 0.5];  % 第四条线，颜色为黑色
            h(2).Color = [0.5, 0.5, 0.5];  % 第五条线，颜色为品红色
            h(1).Color = [0.5, 0.5, 1];  % 第五条线，颜色为品红色
        elseif  curr_day==4||curr_day==5
            h(1).Color = [1, 0, 0];  % 第一条线，颜色为红色
            h(2).Color = [0, 1, 0];  % 第二条线，颜色为绿色
            h(3).Color = [1, 0.5, 0.5];  % 第三条线，颜色为蓝色
            h(4).Color = [0.5, 0.5, 0.5];  % 第四条线，颜色为黑色
            h(5).Color = [0, 0, 0];  % 第五条线，颜色为品红色
            h(6).Color = [1, 0, 0];  % 第五条线，颜色为品红色
            h(7).Color = [0, 1, 0];  % 第五条线，颜色为品红色
            h(8).Color = [1, 0.5, 0.5];  % 第五条线，颜色为品红色
        end
    
    ylim([0 0.4])
    title(['choice day ' num2str(curr_day)])


end

saveas(gcf, fullfile(Path,[ animal '\buffer_image\' animal 'neuron' num2str(curr_cell) 'perievent.jpg']),'jpg')

end








figure
for curr_cell=1:size(data_imgaing_all{curr_day},1)
    imaging1=data_imgaing_all{curr_day}(curr_cell,:)>50;
    peri_event_sample_imaging{curr_cell}=imaging1(sample_event);
    nexttile;
    imagesc(peri_event_sample_imaging{curr_cell}(order,:));
    % plot(mean(peri_event_imging{ii},1))
    % clim([0 1])
    title(num2str(curr_cell))
end





imaging1(sample_event)


figure; hold on
imagesc(track_sample_arm)
plot(round(choice_x1), round(choice_y1))

