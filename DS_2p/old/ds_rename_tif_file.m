% data=dir('G:\CA3_rawdata\CA3_2p\data\1309\data_2p_cell')

% 指定主文件夹路径
mainFolder = 'G:\CA3_rawdata\CA3_2p\data\1675\data_2p_cell'; % 替换为主文件夹路径
% 指定目标文件夹路径
alignFolder = fullfile(mainFolder, 'cell_video_align');
if ~isfolder(alignFolder)
    mkdir(alignFolder);
end
% 获取所有子文件夹
subFolders = dir(mainFolder);
subFolders = subFolders([subFolders.isdir]); % 筛选出文件夹

% 按日期格式筛选文件夹
folderNames = {subFolders.name};
dateFolders = folderNames(~cellfun('isempty', regexp(folderNames, '^\d{4}-\d{2}-\d{2}', 'once')));

% % 按日期顺序排序
% [~, sortIdx] = sort(datetime(dateFolders, 'InputFormat', 'yyyy-MM-dd'));
% dateFolders = dateFolders(sortIdx);

% 遍历每个子文件夹
for i = 1:length(dateFolders)
    currentFolder = fullfile(mainFolder, dateFolders{i}, 'CellVideo2', 'ProcessedVideo');
    
    if isfolder(currentFolder)
        % 获取当前文件夹中的 tif 文件
        tifFiles = dir(fullfile(currentFolder, '*.tif'));
        
        for j = 1:length(tifFiles)
            oldName = fullfile(currentFolder, tifFiles(j).name);
            newName = fullfile(currentFolder, sprintf('%d%s', i, tifFiles(j).name));
            
            % 重命名文件
            movefile(oldName, newName);
            copyfile(newName, alignFolder);
        end
    end
end

disp('重命名完成！');
