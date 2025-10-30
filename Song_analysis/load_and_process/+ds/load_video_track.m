


load(fullfile(ds.locations.filename('server',animal,rec_day,'video_track') ,'mice_position_re.mat'));



% 计算速度
position_speed = nan(size(position_re_X,1)-1,1);
main_preload_vars = who;

dx = diff(position_re_X);
dy = diff(position_re_Y);
dt = diff(position_timelite);
valid = all(~isnan([dx dy dt]), 2);
position_speed(valid) = sqrt(dx(valid).^2 + dy(valid).^2) ./ dt(valid);
position_speed = [position_speed; position_speed(end)];
clearvars('-except',main_preload_vars{:});
