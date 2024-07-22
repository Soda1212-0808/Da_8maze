clear all
% 定义包含CSV文件的文件夹路径
folderPath = 'G:\BaiduNetdiskDownload\data_2p_1\data_path_DLC';
% 获取文件夹中所有以filtered结尾的CSV文件
csvFiles = dir(fullfile(folderPath, '*filtered.csv'))';
mp4Files=dir(fullfile(folderPath, '*.mp4'));
% 获取所有文件的完整路径
csvFilePaths = fullfile({csvFiles.folder}, {csvFiles.name});
mp4FilePaths = fullfile({mp4Files.folder}, {mp4Files.name});

% 使用 cellfun 和 readtable 读取所有CSV文件
tables = cellfun(@readtable, csvFilePaths, 'UniformOutput', false)';
% 使用 vertcat 将所有表格竖向连接起来
mergedTable = vertcat(tables{:});
% 保存合并后的表格到新的CSV文件
writetable(mergedTable, fullfile(folderPath, 'mergedFile.csv'));

X=table2array(mergedTable(:,2))+100;
Y=table2array(mergedTable(:,3))+100;


% 创建一个 VideoReader 对象
v = VideoReader(mp4FilePaths{1});

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

% 
% validIdx=~isnan(X_interp);
% interpIdx = 1:length(X);
% % 插值X坐标
% if any(~validIdx)
%     X_interp(~validIdx) = interp1(interpIdx(validIdx), X_interp(validIdx), interpIdx(~validIdx), 'linear');
%     Y_interp(~validIdx) = interp1(interpIdx(validIdx), Y_interp(validIdx), interpIdx(~validIdx), 'linear');
% end
% 
% figure;
%  imagesc(mask); hold on
% scatter(X_interp,Y_interp)
% plot(X_interp,Y_interp)
