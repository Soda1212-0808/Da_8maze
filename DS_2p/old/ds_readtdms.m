folderA = 'G:\BaiduNetdiskDownload\data_2p_1\data_2p_cell';
cellFiles = dir(fullfile(folderA, '**', 'CellVideo2', 'CellVideo_CHB_Info.tdms'));
miceFiles = dir(fullfile(folderA, '**', 'MiceVideo1', 'MiceVideo_Info.tdms'));
% 提取文件路径
cell_filePaths = fullfile({cellFiles.folder}, {cellFiles.name});
mice_filePaths = fullfile({miceFiles.folder}, {miceFiles.name});
% 定义一个匿名函数，用于读取每个TDMS文件
readTdms = @(file) struct('data', {tdmsread(file)}, 'filePath', file);
% 使用 arrayfun 来读取所有文件数据
cellData(:,1) = arrayfun(readTdms, cell_filePaths, 'UniformOutput', false)';
cellData(:,2) = arrayfun(readTdms, mice_filePaths, 'UniformOutput', false)';
aa=readTdms(mice_filePaths{1})