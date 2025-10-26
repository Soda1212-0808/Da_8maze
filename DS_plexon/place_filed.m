%% 参数设置
bin_size = 3;  % 位置 bin 大小
occupancy_smoothing = 1;  % occupancy map 高斯平滑核宽度
shuffle_num = 1000;  % shuffle 次数
min_field_area = 9;  % 连通区域中，位置野的最小像素数
alpha = 0.05;  % 显著性水平

%% Step 1. 构建位置 bin
% x_edges = min(position(:,1)) : bin_size : max(position(:,1));
% y_edges = min(position(:,2)) : bin_size : max(position(:,2));
[xx, yy] = meshgrid(x_edges(1:end-1)+bin_size/2, y_edges(1:end-1)+bin_size/2);

%% Step 2. spike 对齐位置
spike_frame_idx = knnsearch(position_timelite, spike_times);
spike_pos = position(spike_frame_idx, :);

%% Step 3. occupancy map 和 spike map
occupancy_map = histcounts2(position(:,1), position(:,2), x_edges, y_edges);
spike_map = histcounts2(spike_pos(:,1), spike_pos(:,2), x_edges, y_edges);

%% Step 4. 平滑处理
g = fspecial('gaussian', [3 3], occupancy_smoothing);
smoothed_occupancy = conv2(occupancy_map, g, 'same');
smoothed_spike = conv2(spike_map, g, 'same');

%% Step 5. ratemap
ratemap = smoothed_spike ./ (smoothed_occupancy + eps);

%% Step 6. 空间信息量 SI 计算
p_ij = smoothed_occupancy / sum(smoothed_occupancy(:));
lambda_ij = ratemap;
lambda_bar = sum(p_ij(:) .* lambda_ij(:));
info_map = p_ij .* (lambda_ij ./ lambda_bar) .* log2(lambda_ij ./ lambda_bar + eps);
SI_real = nansum(info_map(:));

%% Step 7. shuffle 检验 SI
SI_shuffled = zeros(shuffle_num,1);
total_time = position_time(end) - position_time(1);
for i = 1:shuffle_num
    shift = rand * total_time;
    spike_times_shuffled = mod(spike_times + shift, total_time);
    shuffled_idx = knnsearch(position_time, spike_times_shuffled);
    shuffled_pos = position(shuffled_idx, :);
    shuffled_spike_map = histcounts2(shuffled_pos(:,1), shuffled_pos(:,2), x_edges, y_edges);
    shuffled_spike_smooth = conv2(shuffled_spike_map, g, 'same');
    shuffled_ratemap = shuffled_spike_smooth ./ (smoothed_occupancy + eps);
    
    lambda_ij = shuffled_ratemap;
    lambda_bar = sum(p_ij(:) .* lambda_ij(:));
    info_map = p_ij .* (lambda_ij ./ lambda_bar) .* log2(lambda_ij ./ lambda_bar + eps);
    SI_shuffled(i) = nansum(info_map(:));
end

% 判断是否为位置细胞
p_value_SI = mean(SI_shuffled >= SI_real);
is_place_cell = p_value_SI < alpha;

%% Step 8. shuffle 检验位置野阈值
max_shuffled = zeros(shuffle_num,1);
for i = 1:shuffle_num
    shift = rand * total_time;
    shuffled_times = mod(spike_times + shift, total_time);
    shuffled_idx = knnsearch(position_time, shuffled_times);
    shuffled_pos = position(shuffled_idx, :);
    shuffled_spike_map = histcounts2(shuffled_pos(:,1), shuffled_pos(:,2), x_edges, y_edges);
    shuffled_spike_smooth = conv2(shuffled_spike_map, g, 'same');
    shuffled_ratemap = shuffled_spike_smooth ./ (smoothed_occupancy + eps);
    max_shuffled(i) = max(shuffled_ratemap(:));
end
firing_threshold = prctile(max_shuffled, 95);

%% Step 9. 提取位置野
has_place_field = false;
place_field_mask = false(size(ratemap));
if is_place_cell
    binary_map = ratemap > firing_threshold;
    CC = bwconncomp(binary_map, 8);
    for k = 1:CC.NumObjects
        if numel(CC.PixelIdxList{k}) >= min_field_area
            place_field_mask(CC.PixelIdxList{k}) = true;
        end
    end
    has_place_field = any(place_field_mask(:));
end

%% Step 10. 输出 & 可视化
fprintf('Spatial Information: %.3f bits/spike\n', SI_real);
fprintf('p-value (SI shuffle): %.4f\n', p_value_SI);
fprintf('位置细胞: %s\n', ternary(is_place_cell, '✅ 是', '❌ 否'));
fprintf('位置野: %s\n', ternary(has_place_field, '✅ 有显著位置野', '❌ 无'));

figure;
imagesc(xx(1,:), yy(:,1), ratemap'); axis xy image;
colormap hot; colorbar;
title(sprintf('Ratemap (SI = %.2f, p = %.4f)', SI_real, p_value_SI));
hold on;
if has_place_field
    contour(xx, yy, place_field_mask', [0.5 0.5], 'LineWidth', 2, 'Color', 'r');
end

%% 工具函数 ternary
function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
