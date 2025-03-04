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
days = cellfun(@(x) x(1:10), fileNames, 'UniformOutput', false)';

animal_path = cellfun(@readtable, csvFilePaths, 'UniformOutput', false)';

save(fullfile(Path,animal,newfolderName,['merged_mice_path.mat']),'animal_path','days')

end