   
    %获取MP4文件
    mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.mp4'));
    if isempty(mp4Files)
        mp4Files=dir(fullfile(Path, animal ,'data_path_DLC' ,'*.avi'));
    end
    mp4FilePaths=fullfile({mp4Files.folder}, {mp4Files.name});
    % 创建一个 VideoReader 对象
    v = VideoReader(mp4FilePaths{1});
    framerate=v.framerate;
    % 读取第一帧
    firstFrame = readFrame(v);
    [rows, cols, channels] = size(firstFrame);
    % 创建一个新图像，尺寸比原始图像大200像素（上下左右各100像素）
    newRows = rows + 200;
    newCols = cols + 200;
    paddedFrame = uint8(255 * ones(newRows, newCols, channels)); % 白色背景
    % 将原始图像复制到新图像的中心
    paddedFrame(101:100+rows, 101:100+cols, :) = firstFrame;

    % 若之前未绘制目标区域并保存文件，在轨迹图像上绘制目标区域
    if exist(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'))~=2
        display_next_frame_on_scroll(mp4FilePaths{1})
        figure;
        % 显示填充后的第一帧
        imshow(paddedFrame);
        hold on
        cellfun(@(x,y) scatter(x,y),X_position,Y_position,'UniformOutput',false)

        % plot(X,Y)
        % mask=roipoly;
        % 指定需要绘制的多边形区域数量
        numPolygons = 1; % 你可以根据需要改变此值
        BW1=cell(1,numPolygons );
        for k = 1:numPolygons
            % 绘制多边形区域
            BW = roipoly;
            BW1{k}=BW;
            % 将当前多边形区域添加到组合掩码中
            %     BW_combined = BW_combined | BW;

            % 显示当前多边形区域的边界
            boundary = bwboundaries(BW);
            plot(boundary{1}(:,2), boundary{1}(:,1), 'LineWidth', 2);
            mark_x = boundary{1}(1, 2);
            mark_y = boundary{1}(1, 1);
            text(mark_x, mark_y, num2str(k-1), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');

        end
        saveas(gcf,fullfile(Path, animal ,bufferfolderName,'grab_picture.jpg'),'jpeg')
        save(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'),'BW1','-mat')
        % save(fullfile(Path, animal , bufferfolderName,'grab_picture.mat'),'BW1','recordedFrameCount','-mat')
    else
        load(fullfile(Path, animal ,bufferfolderName,'grab_picture.mat'))
    end
