clear all
% 定义包含CSV文件的文件夹路径
Path = 'D:\SDdata\data_2p_1';

 Path = 'G:\CA3_rawdata\CA3_2p\data';
 animals={'1306','1307','1309','1311','1312','1646','1974','1976'};
 
for curr_animal=1:8
    animal=animals{curr_animal};

DLC_folder='data_path_DLC';
newfolderName = 'merged file';
if exist(fullfile(Path,newfolderName), 'dir') ~= 7
    mkdir(fullfile(Path,newfolderName));
    disp(['Folder "', newfolderName, '" created.']);
end
% 获取文件夹中所有以filtered结尾的CSV文件
csvFiles = dir(fullfile(Path,animal,DLC_folder, '*filtered.csv'))';

% 获取所有CSV文件的完整路径
csvFilePaths = fullfile({csvFiles.folder}, {csvFiles.name});

% 提取文件名
fileNames = {csvFiles.name};

% 提取文件名的第5到第14位字符
groupIDs = cellfun(@(x) x(1:10), fileNames, 'UniformOutput', false);

% 找到唯一的组ID
[uniqueGroupIDs, ~, groupIdx] = unique(groupIDs);

% 按组读取并合并文件
cell_animal_path = cell([size(uniqueGroupIDs,2),1]); % 预分配单元数组
animal_path=struct;
for id = 1:numel(uniqueGroupIDs)
    % 读取当前组的所有文件
    tables = cellfun(@readtable, csvFilePaths(groupIdx == id), 'UniformOutput', false);

    % 获取第一个 table 的列名作为标准
    standardColumnNames = tables{1}.Properties.VariableNames;

% 重命名所有 table 的列名以匹配标准列名
tables = cellfun(@(tbl) setfield(tbl, 'Properties', ...
    setfield(tbl.Properties, 'VariableNames', standardColumnNames)), ...
    tables, 'UniformOutput', false);
    % 合并表格
    cell_animal_path{id} = vertcat(tables{:});

%     animal_path.(['date_' strrep(uniqueGroupIDs{id}, '-', '_')])=vertcat(tables{:});
end


save(fullfile(Path,newfolderName,['merged_mice_path.mat']),'cell_animal_path','animal_path')

% % 保存每个组的合并表格到新的CSV文件
% cellfun(@(table, id) ...
%     writetable(table, fullfile(fullfile(folderPath,newfolderName), ['mergedGroup_' id '.csv'])), ...
%     animal_path, uniqueGroupIDs, 'UniformOutput', false);
end