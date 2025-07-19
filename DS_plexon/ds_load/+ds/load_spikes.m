% load spikes
neuron_files=dir(fullfile(Path, animal , rec_day,'recording' ,'*.t'));
spikes_all = arrayfun(@(f) readmclusttfile(fullfile(f.folder, f.name))'/10000, ...
    neuron_files, 'UniformOutput', false);
spike_freq=cell2mat(cellfun(@(x) length(x)/x(end),spikes_all,'UniformOutput',false));

templates_all = cell2mat(arrayfun(@(idx, n) repmat(idx, n, 1),...
    (1:numel(cellfun(@numel, spikes_all)))', cellfun(@numel, spikes_all), 'UniformOutput', false));
[spike_times_timelite, sort_idx] = sort(cell2mat(spikes_all));
spike_templates = templates_all(sort_idx);
