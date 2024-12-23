clear all
% 设置主文件夹路径
Path = 'E:\data_8_maze\data_2p_1\1464';
data_2p_folder='data_2p_cell';
newfolderName = 'merged file';
if exist(fullfile(Path,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end
% 获取主文件夹中的所有子文件夹
subfolders = dir(fullfile(Path,data_2p_folder));
subfolderNames = {subfolders.name};

% 过滤出日期格式的子文件夹
datePattern = '\d{4}-\d{2}-\d{2}';
dateFolders = subfolderNames(cellfun(@(x) ~isempty(regexp(x, datePattern, 'once')) && isfolder(fullfile(Path,data_2p_folder, x, 'MiceVideo1')), subfolderNames));
% 获取所有日期
dates = cellfun(@(x) regexp(x, datePattern, 'match', 'once'), dateFolders, 'UniformOutput', false);

dates_maze = cellfun(@(x) regexp(x, '8-maze', 'match', 'once'), dateFolders, 'UniformOutput', false);

for i=1:length(dates_maze)

load(fullfile(Path,data_2p_folder, dateFolders{i}, 'matchTable.mat'))
if ~isempty(dates_maze{i})
 
matchTable =[matchTable num2cell(ones(size(matchTable,1), 1))];
else matchTable =[matchTable num2cell(2*ones(size(matchTable,1), 1))];
    
end
save(fullfile(Path,data_2p_folder, dateFolders{i}, 'matchTable_new.mat'),'matchTable')
end





% 去重日期
uniqueDates = unique(dates);

% 创建一个结构体数组，用于按日期存储tdms数据
animal_timepoint = struct();
cell_timepoint = struct();


% 读取并合并同一日期的tdms文件
for dateIdx = 1:numel(uniqueDates)
    dateStr = uniqueDates{dateIdx};
    sameDateFolders = dateFolders(strcmp(dates, dateStr));
    
    %读取并合并同一天的matchTable文件
    match_table_mice = arrayfun(@(x) load(fullfile(Path,data_2p_folder, x{1}, 'matchTable_new.mat')), sameDateFolders, 'UniformOutput', false);
     % 合并数据
    buffer3 = vertcat(match_table_mice{:});
    tablesArray = cellfun(@(x) x.matchTable, num2cell(buffer3), 'UniformOutput', false);
    combined_match_table_mice = vertcat(tablesArray{:});
    % 存储到结构体
    animal_match_table.(['date_' strrep(dateStr, '-', '_')]) = combined_match_table_mice;


    % 读取并合并同一天的tdms文件
    tdmsDataList_mice = arrayfun(@(x) tdmsread(fullfile(Path,data_2p_folder, x{1}, 'MiceVideo1', 'MiceVideo_Info.tdms')), sameDateFolders, 'UniformOutput', false);
    % 合并数据
    buffer1 = vertcat(tdmsDataList_mice{:});
    combinedTdmsData_mice = vertcat(buffer1{:});

    combinedTdmsData_mice.timepoint = seconds(datetime(combinedTdmsData_mice.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime(combinedTdmsData_mice.Time(1), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS'));
    combinedTdmsData_mice.ID=str2double(combinedTdmsData_mice.ID);
    % 存储到结构体
    animal_timepoint.(['date_' strrep(dateStr, '-', '_')]) = combinedTdmsData_mice;

    tdmsDataList_cell = arrayfun(@(x) tdmsread(fullfile(Path,data_2p_folder, x{1}, 'CellVideo2', 'CellVideo_CHB_Info.tdms')), sameDateFolders, 'UniformOutput', false);
    buffer2= vertcat(tdmsDataList_cell{:});
    combinedTdmsData_cell = vertcat(buffer2{:,2});
    combinedTdmsData_cell.Time = regexprep(combinedTdmsData_cell.Time, '-(?=[^-]*$)', ':');
    combinedTdmsData_cell.Time = regexprep(combinedTdmsData_cell.Time, '-(?=[^-]*$)', ':');

    % 将字符串转换为 datetime 类型
    combinedTdmsData_cell.timepoint = seconds(datetime(combinedTdmsData_cell.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS')-datetime(combinedTdmsData_cell.Time(1), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS'));
    combinedTdmsData_cell.ID=str2double(combinedTdmsData_cell.ID);

    cell_timepoint.(['date_' strrep(dateStr, '-', '_')]) = combinedTdmsData_cell;


end

save(fullfile(Path,newfolderName,['merged_mice_cell_timepoint.mat']),'animal_timepoint','cell_timepoint','animal_match_table')

