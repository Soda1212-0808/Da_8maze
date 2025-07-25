clear all
% 设置主文件夹路径
%  Path = 'E:\data_8_maze\data_2p_1\1464';
Path = 'G:\CA3_rawdata\CA3_2p\data';
animals={'1306','1307','1309','1311','1312','1646','1974','1976'};

for curr_animal=1:length(animals)
    preload_vars = who;

    animal=animals{curr_animal};
    data_2p_folder='data_2p_cell';
    newfolderName = 'merged file';
    if exist(fullfile(Path,animal,newfolderName), 'dir') ~= 7
        mkdir(fullfile(Path,animal,newfolderName));
        disp(['Folder "', newfolderName, '" created.']);
    end
    % 获取主文件夹中的所有子文件夹
    subfolders = dir(fullfile(Path,animal,data_2p_folder));
    subfolderNames = {subfolders.name};

    % 过滤出日期格式的子文件夹
    datePattern = '\d{4}-\d{2}-\d{2}';
    dateFolders = subfolderNames(cellfun(@(x) ~isempty(regexp(x, datePattern, 'once')) && isfolder(fullfile(Path,animal,data_2p_folder, x, 'MiceVideo1')), subfolderNames));
    % 获取所有日期
    dates = cellfun(@(x) regexp(x, datePattern, 'match', 'once'), dateFolders, 'UniformOutput', false);

    dates_maze = cellfun(@(x) regexp(x, '8-maze', 'match', 'once'), dateFolders, 'UniformOutput', false);

    for i=1:length(dates_maze)

        load(fullfile(Path,animal,data_2p_folder, dateFolders{i}, 'matchTable.mat'))
        if ~isempty(dates_maze{i})

            matchTable =[matchTable num2cell(ones(size(matchTable,1), 1))];
        else matchTable =[matchTable num2cell(2*ones(size(matchTable,1), 1))];

        end
        save(fullfile(Path,animal,data_2p_folder, dateFolders{i}, 'matchTable_new.mat'),'matchTable')
    end

    % 去重日期
    all_Dates = unique(dates)';

    % 创建一个结构体数组，用于按日期存储tdms数据

    %     animal_match=cell(numel(all_Dates),1);
    %     animal_timepoint =cell(numel(all_Dates),1);
    %     cell_timepoint = cell(numel(all_Dates),1);

    % 读取并合并同一日期的tdms文件
    for dateIdx = 1:numel(all_Dates)
        dateStr = all_Dates{dateIdx};

        sameDateFolders = dateFolders(strcmp(dates, dateStr));

        %读取并合并同一天的matchTable文件
        match_table_mice = arrayfun(@(x) load(fullfile(Path,animal,data_2p_folder, x{1}, 'matchTable_new.mat')), sameDateFolders, 'UniformOutput', false);
        % 合并数据
        buffer3 = vertcat(match_table_mice{:});
        tablesArray = cellfun(@(x) x.matchTable, num2cell(buffer3), 'UniformOutput', false);
        combined_match_table_mice = vertcat(tablesArray{:});
        % 存储到结构体

        %         animal_match{dateIdx}=combined_match_table_mice;
        animal_match=combined_match_table_mice;

        % 读取并合并同一天的tdms文件
        tdmsDataList_mice = arrayfun(@(x) tdmsread(fullfile(Path,animal,data_2p_folder, x{1}, 'MiceVideo1', 'MiceVideo_Info.tdms')), sameDateFolders, 'UniformOutput', false);
        % 合并数据
        buffer1 = vertcat(tdmsDataList_mice{:});
        combinedTdmsData_mice = vertcat(buffer1{:});

        combinedTdmsData_mice.timepoint = seconds(datetime(combinedTdmsData_mice.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime(combinedTdmsData_mice.Time(1), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS'));
        combinedTdmsData_mice.ID=str2double(combinedTdmsData_mice.ID);
        % 存储到结构体

        %         animal_timepoint{dateIdx}= combinedTdmsData_mice;
        video_timepoint= combinedTdmsData_mice;

        tdmsDataList_cell = arrayfun(@(x) tdmsread(fullfile(Path,animal,data_2p_folder, x{1}, 'CellVideo2', 'CellVideo_CHB_Info.tdms')), sameDateFolders, 'UniformOutput', false);
        buffer2= vertcat(tdmsDataList_cell{:});
        combinedTdmsData_cell = vertcat(buffer2{:,2});
        combinedTdmsData_cell.Time = regexprep(combinedTdmsData_cell.Time, '-(?=[^-]*$)', ':');
        combinedTdmsData_cell.Time = regexprep(combinedTdmsData_cell.Time, '-(?=[^-]*$)', ':');

        % 将字符串转换为 datetime 类型
        combinedTdmsData_cell.timepoint = seconds(datetime(combinedTdmsData_cell.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime(combinedTdmsData_cell.Time(1), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS'));
        combinedTdmsData_cell.ID=str2double(combinedTdmsData_cell.ID);

%         cell_timepoint{dateIdx}=combinedTdmsData_cell;
        cell_timepoint=combinedTdmsData_cell;
        if exist(fullfile(Path,animal,dateStr), 'dir') ~= 7
        mkdir(fullfile(Path,animal,dateStr));
        disp(['Folder "', dateStr, '" created.']);
    end
   save(fullfile(Path,animal,dateStr,[dateStr ,'_video_cell_match.mat']),'animal_match','video_timepoint','cell_timepoint','-v7.3')

    end

    % save(fullfile(Path,animal,newfolderName,['merged_mice_cell_timepoint.mat']),'uniqueDates','animal_timepoint','cell_timepoint','animal_match_table','animal_match_cell','animal_timepoint_cell','cell_timepoint_cell')
%     save(fullfile(Path,animal,newfolderName,['merged_mice_cell_timepoint.mat']),'all_Dates','animal_match','animal_timepoint','cell_timepoint')


    clearvars('-except',preload_vars{:});

end