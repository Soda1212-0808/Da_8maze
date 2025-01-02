

track=BW1{2};


for i=1:size(all_data_event.all_event_frame{1},1)
    range=all_data_event.all_event_frame{1}(i,1:2);
    x1=X_filter_speed_all{1, 1}(range(1):range(2));
    y1=Y_filter_speed_all{1, 1}(range(1):range(2));
    idx=~isnan(x1);
    in_target_area=nan(length(x1),1);
    in_target_area(idx) = track(sub2ind(size(track), round(y1(idx)), round(x1(idx))));
    if ~isempty(find(in_target_area==1,1,'last'))
        result(i)=range(1)-1+ find(in_target_area==1,1,'last');
    else  result(i)=range(1)
    end
end







figure; hold on
imagesc(track)
plot(round(x1(idx)), round(y1(idx)))

