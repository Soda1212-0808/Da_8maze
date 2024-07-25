clc
clear;
% delete(gcp('nocreate'));
% numCore = feature('numcores');
% parpool(numCore-1);
pre=-5;
post=5;
width=0.05;
bin=0.05;
fr=1/bin;
step=(post-pre)/bin;

Path = 'H:\CA3_reprocess\each mice\';    % 设置数据存放的文件夹路径
animals={'DCA3-9','DCA3-10','DCA3-11','DCA3-12','DCA3-14','DCA3-17','DCA3-20'};
for curr_animal=1:length(animals)
    % 获取文件夹中的所有内容
    contents = dir([Path animals{curr_animal}]);
    % 获取所有子文件夹的名称
    recording_files = {contents(([contents.isdir] & ~ismember({contents.name}, {'.', '..'}))).name};

    for curr_file=1:length(recording_files)

        nex5_file=dir(fullfile([Path animals{curr_animal} '\' recording_files{curr_file} ],'*.nex5'));
        file_name= {nex5_file.name}';
        mkdir([Path animals{curr_animal} '\' recording_files{curr_file} '\spikes\']);
        trial=struct;
        data=readNex5File([nex5_file.folder '\' nex5_file.name]);
        neuron_files=dir(fullfile([Path animals{curr_animal} '\' recording_files{curr_file} ],'*.t64'));
        events=length(data.events);
        n_c_b=struct;
        dat=struct;
        for current_neuron=1:length(neuron_files)

            neuron_name=neuron_files(current_neuron).name(1:end-4);
            spike_whole=readmclusttfile( [neuron_files(current_neuron).folder '\' neuron_files(current_neuron).name])'/10000;
            c_b=struct;

            for current_event=1:events
                t=data.events{current_event};

                c_b(current_event).name=t.name;
                c_b(current_event).count_bin=[];

                t_s=t.timestamps;
                p1=~(strcmp(c_b(current_event).name,'Start') || strcmp(c_b(current_event).name,'Stop'));
                p2=~isempty(t_s);
                if p1&p2

                    % 计算每个事件的尖峰位置
                    num_events = length(t_s);
                    spike_positions = arrayfun(@(t) spike_whole(spike_whole > (t - 5) & spike_whole < (t + 5)) - t, t_s, 'UniformOutput', false);
                    % 存储尖峰位置
                    c_b(current_event).spike_position = spike_positions;

                    % 计算每个 bin 内的尖峰数量
                    % 构建 bin 的边界
                    bin_edges = repmat(t_s, 1, step + 1) + pre + bin * (0:step) + [zeros(num_events, 1), repmat(width, num_events, step)];
                    % 将边界转换为一个单元数组
                    edges_cell = num2cell(bin_edges, 2);
                    % 计算每个 bin 内的尖峰数量
                    counts = cellfun(@(e) histcounts(spike_whole, e), edges_cell, 'UniformOutput', false);
                    % 将结果转换为矩阵
                    c_b(current_event).count_bin = cell2mat(counts)';

                    c_b(current_event).firing_rate= c_b(current_event).count_bin/width;
                    c_b(current_event).firing_rate_gaussian=smoothdata(c_b(current_event).firing_rate,1, "gaussian", 11);
                    c_b(current_event).firing_rate_mean= mean(c_b(current_event).firing_rate,2);
                    c_b(current_event).firing_rate_gaussian_mean=smoothdata(c_b(current_event).firing_rate_mean, "gaussian", 11);


                end
            end
            %         save([Path newfile name_p(1:end-5) '_' neuron_name '.mat'],'c_b', '-v7.3');
            save([Path animals{curr_animal} '\' recording_files{curr_file} '\spikes\'   neuron_name '.mat'],'c_b');


        end

    end




end
