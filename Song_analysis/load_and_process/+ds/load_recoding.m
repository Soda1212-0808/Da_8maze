% ap.load_recording
%
% Load and prepare all/selected parts of a recording

warning on

%% Check that the necessary parameters are included

% If no animal - choose from list of animals
if ~exist('animal','var') || isempty(animal)
    if ~exist('rec_day','var') || isempty(rec_day)
        % (if no day - return all animals)
        data_dir = dir(plab.locations.local_data_path);
        animals_all = {data_dir(3:end).name};
    else
        % (if day specified - return animals with those days)
        day_dir = dir(fullfile(plab.locations.server_data_path,'*',rec_day));
        animals_all = erase(unique({day_dir.folder})',{plab.locations.server_data_path,rec_day,filesep});
    end
    animal_idx = listdlg('PromptString','Select animal:', ...
        'ListString',animals_all,'ListSize',[300,200], ...
        'SelectionMode','single');
    animal = animals_all{animal_idx};
    verbose = true;
end

% If no day - choose from list of days
if ~exist('rec_day','var') || isempty(rec_day)
    animal_recordings = ds.find_recordings(animal);
    recording_days = {animal_recordings.day};
    day_idx = listdlg('PromptString','Select day:', ...
        'ListString',recording_days,'ListSize',[300,200], ...
        'SelectionMode','single');
    rec_day = recording_days{day_idx};
    verbose = true;
end

% % If no time - choose from list of workflows
% if ~exist('rec_time','var') || isempty(rec_time)
%     recordings = plab.find_recordings(animal,rec_day);
%     % (use only recordings with Bonsai workflows)
%     bonsai_recordings = find(~cellfun(@isempty,[recordings.workflow]));
%     rec_idx = listdlg('PromptString','Select workflow:', ...
%         'ListString',recordings.workflow(bonsai_recordings),'ListSize',[300,200], ...
%         'SelectionMode','single');
%     rec_time = recordings.recording{bonsai_recordings(rec_idx)};
%     verbose = true;
% end


%% Define what to load

if ~exist('verbose','var')
    verbose = false;
end

if verbose; fprintf('Loading %s, %s, Recording %s\n', animal, rec_day); end;

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts','var')
    load_parts.video_track = true;
    load_parts.ca_2p = true;
    load_parts.tetrode = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'video_track')
        load_parts.video_track = false;
    end
    if ~isfield(load_parts,'ca_2p')
        load_parts.ca_2p = false;
    end
    if ~isfield(load_parts,'tetrode')
        load_parts.tetrode = false;
    end
end



%% Load experiment components

% Load timelite and associated inputs
ds.load_events



% Load mousecam
if load_parts.video_track
    ds.load_video_track
end

% Load 2p
if load_parts.ca_2p && ...
        exist(ds.locations.filename('local',animal,rec_day,'ca_2p_recording'),'dir')
    ap.load_ca_2p
end

% Load tetrode
if load_parts.tetrode && ...
        exist(ds.locations.filename('local',animal,rec_day,'tetrode_recording'),'dir')
    ds.load_spikes
end

if verbose; disp('Finished.'); end;













