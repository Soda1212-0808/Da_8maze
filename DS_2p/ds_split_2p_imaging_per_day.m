clear all
% 定义包含CSV文件的文件夹路径
Path = 'G:\CA3_rawdata\CA3_2p\data';    % 设置数据存放的文件夹路径
% animals={'1464'};
animals={'1646','1306','1307','1309','1311','1312','1974','1976'};
for curr_animal=1:length(animals)
        animal=animals{curr_animal};
        

        data_imaging=load(fullfile(Path,animal,'data_2p_cell\cell_video_align\suite2p\plane0','Fall.mat'));

        [uniqueDates, ~, fie_idx] = unique(regexp(strtrim(string(data_imaging.ops.filelist) )  , '\d{4}-\d{2}-\d{2}', 'match', 'once'), 'stable');
        frames_per_day = accumarray(fie_idx,  data_imaging.ops.frames_per_file');
        valided_cells=data_imaging.spks(data_imaging.iscell(:,1)==1,:);
        valided_cells_per_day = mat2cell(valided_cells, size(valided_cells,1), frames_per_day)';
        for curr_day=1:length(valided_cells_per_day)

        spikes=valided_cells_per_day{curr_day};
        save(fullfile(Path,animal,uniqueDates(curr_day),'image_2p/',strcat(uniqueDates(curr_day) ,'_spikes.mat')),'spikes','-v7.3')
        end
end







