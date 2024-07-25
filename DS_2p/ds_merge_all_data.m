clear all

Path = 'G:\BaiduNetdiskDownload\data_2p_1\merged file';

files=dir(fullfile(Path,'*.mat'))
load(fullfile(files(1).folder,files(1).name))
load(fullfile(files(2).folder,files(2).name))

% 获取animal_match_table的所有字段名（即所有日期）
dateFields = fieldnames(animal_match_table);

% 定义一个匿名函数来对每个日期字段进行操作
operateOnDate = @(dateField) ...
    [animal_match_table.(dateField) ...
    table2cell(animal_path.(dateField)(cell2mat(animal_match_table.(dateField)(:,3)) + 1, :))];

% 使用cellfun对所有日期字段进行操作，并将结果存入newdata结构体中
newdata = cell2struct(cellfun(operateOnDate, dateFields, 'UniformOutput', false), dateFields);

currdata=newdata.(dateFields{1});
X=cell2mat(currdata(:,8))+100;
Y=cell2mat(currdata(:,9))+100;


folderPath = 'G:\BaiduNetdiskDownload\data_2p_1\data_path_DLC';
mp4Files=dir(fullfile(folderPath, '*.mp4'));

currindx=  find(contains({mp4Files.name},  strrep(dateFields{1}(6:end), '_', '-')),1);
mp4FilePaths=fullfile({mp4Files.folder}, {mp4Files.name});
% 创建一个 VideoReader 对象
v = VideoReader(mp4FilePaths{currindx});

% 读取第一帧
firstFrame = readFrame(v);
[rows, cols, channels] = size(firstFrame);
% 创建一个新图像，尺寸比原始图像大200像素（上下左右各100像素）
newRows = rows + 200;
newCols = cols + 200;
paddedFrame = uint8(255 * ones(newRows, newCols, channels)); % 白色背景
% 将原始图像复制到新图像的中心
paddedFrame(101:100+rows, 101:100+cols, :) = firstFrame;


figure
% 显示填充后的第一帧
imshow(paddedFrame);

hold on
scatter(X,Y)
% plot(X,Y)
mask=roipoly;


X_interp=X;
Y_interp=Y;

% 找到 mask 中对应位置为 0 的索引
invalid_indices = mask(sub2ind(size(mask), fix(Y_interp+1), fix(X_interp+1))) == 0;
% 将对应的 X_interp 和 Y_interp 位置设为 NaN
X_interp(invalid_indices) = NaN;
Y_interp(invalid_indices) = NaN;


figure;
nexttile;
 imagesc(mask); hold on
scatter(X,Y)
nexttile;
 imagesc(mask); hold on
scatter(X_interp,Y_interp)
plot(X_interp,Y_interp)
