% load events
%处理event数据


event_file=dir(fullfile(ds.locations.filename('server',animal,rec_day,'behavior') ,'*.csv'));

data_event=csvread(fullfile(event_file.folder , event_file.name));

events_name= unique (data_event(:,1));

event_times=cell(8,1);
event_times(1:4) = arrayfun(@(event) arrayfun(@(id) sort(data_event(data_event(:,1) == id, event)),...
    events_name, 'UniformOutput', false),...
    3:6, 'UniformOutput', false);

events_name_swapped = mod(events_name,10)*10 + floor(events_name/10);
event_times(5:8) = arrayfun(@(event) arrayfun(@(id) sort(data_event(data_event(:,1) == id, event)),...
    events_name_swapped, 'UniformOutput', false),...
    7:10, 'UniformOutput', false);
event_labels={'sample begin';'sample arm';'sample reward';'delay';...
    'choice begin';'choice arm';'choice reward';'end'};
events=struct('event_labels',event_labels,'event_times',event_times);

arm_times=arrayfun(@(event) arrayfun(@(id1,id2) sort([ data_event(data_event(:,1)==id1,event);...
    data_event(data_event(:,1)==id2&data_event(:,2)==1,event+4);...
    data_event(data_event(:,1)==id1&data_event(:,2)==0,event+4)]),...
    events_name,events_name_swapped, 'UniformOutput', false),...
    3:6, 'UniformOutput', false);



% 如果是ca_2p,则event是帧模式

allInts = all(mod(data_event(:),1) == 0);
if allInts
    event_file=dir(fullfile(ds.locations.filename('server',animal,rec_day) ,'*video_cell_match.mat'));
    load(fullfile(event_file.folder,event_file.name))




    A=cell2mat(animal_match(:,3))';
    B=str2double(animal_match(:,1))';

    idx = find(diff(A) < 0) + 1; % 找到下降点的索引

    if isempty(idx)
        animal_match(:,8) = num2cell(A');
        animal_match(:,7) = num2cell(B');

    else
        % 计算累积调整量，排除初始的0
        adjustments_A = cumsum([0, A(idx-1)]);

        % 利用广播机制创建调整矩阵，并计算总调整量
        delta_A = max((1:numel(A) >= idx(:)) .* adjustments_A(2:end)', [],1);
        % 应用调整
        animal_match(:,8) = num2cell((A + delta_A)');

        adjustments_B = cumsum([0, B(idx-1)]);
        % 利用广播机制创建调整矩阵，并计算总调整量
        delta_B = max((1:numel(B) >= idx(:)) .* adjustments_B(2:end)',[], 1);
        % 应用调整
        animal_match(:,7) = num2cell((B + delta_B)');
    end

    [~,idxx]=  unique(double(cell2mat(animal_match(:, 8))), 'first') ;
    data_event_to_2p =[data_event(:,1:2) ...
        round(interp1(double(cell2mat(animal_match(idxx, 8))), cell2mat(animal_match(idxx, 7)),data_event(:,3:10), 'linear', 'extrap')) ];

end
