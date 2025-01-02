% 主函数：根据日期合并文件数据
function dataByDate = mergeExcelFilesByDate(folderPath, startCol, endCol, rowStart)
    % 设置文件夹路径
    if nargin < 1
        folderPath = 'C:\Your\Folder\Path'; % 替换为实际路径
    end
    
    % 默认列范围
    if nargin < 2, startCol = 'J'; end
    if nargin < 3, endCol = 'O'; end
    if nargin < 4, rowStart = 2; end

    % 获取文件夹中所有 .xlsx 文件的完整路径
    filePattern = fullfile(folderPath, '*.xlsx');
    xlsxFiles = dir(filePattern);
    filePaths = fullfile({xlsxFiles.folder}, {xlsxFiles.name});
    fileNames = {xlsxFiles.name};

    % 提取日期并分组文件
    fileDates = extractDatesFromFilenames(fileNames);
    uniqueDates = unique(fileDates);

    % 合并数据
    dataByDate = cell2struct(arrayfun(@(date) mergeFilesForDate(filePaths, fileDates, date{1}, startCol, endCol, rowStart), ...
        uniqueDates, 'UniformOutput', false), uniqueDates, 1);
end

% --- 子函数 1：提取日期 ---
function fileDates = extractDatesFromFilenames(fileNames)
    datePattern = '\d{4}_\d{1,2}_\d{1,2}'; % 匹配日期模式
    fileDates = cellfun(@(name) regexp(name, datePattern, 'match', 'once'), fileNames, 'UniformOutput', false);
end

% --- 子函数 2：合并同日期文件数据 ---
function mergedData = mergeFilesForDate(filePaths, fileDates, targetDate, startCol, endCol, rowStart)
    % 筛选出目标日期的文件
    targetFiles = filePaths(strcmp(fileDates, targetDate));

    % 合并数据
    mergedData = vertcat(cellfun(@(file) table2array(readSpecifiedRange(file, startCol, endCol, rowStart)), ...
        targetFiles, 'UniformOutput', false){:});
end

% --- 子函数 3：读取指定范围的数据 ---
function data = readSpecifiedRange(filePath, startCol, endCol, rowStart)
    % 动态生成范围
    rangeStr = dynamicRange(filePath, startCol, endCol, rowStart);
    
    % 读取指定范围的数据
    data = readtable(filePath, 'Range', rangeStr);
end

% --- 子函数 4：动态生成范围 ---
function rangeStr = dynamicRange(filePath, startCol, endCol, rowStart)
    % 获取文件的行数
    [~, ~, raw] = xlsread(filePath);
    maxRow = size(raw, 1);
    
    % 构造范围字符串
    rangeStr = sprintf('%s%d:%s%d', startCol, rowStart, endCol, maxRow);
end
