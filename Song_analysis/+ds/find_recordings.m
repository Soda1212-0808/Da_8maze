function recordings = find_recordings(animal,recording_day,workflow)
% recordings = find_recordings(animal,recording_day,workflow)
%
% Find recordings for a given animal/day/workflow, package info by day.
%
% Input:
% recording_day - (as 'yyyy-mm-yy'), can be multiple in a cell array
% workflow - can be multiple (e.g. {'taskA','taskB'}), can include * as
% wildcard (e.g. 'task*' will return both 'taskA' and 'taskB').
%
%  - if only animal defined: find all recordings across all days
%  - if animal & day defined: find all recordings within one day
%  - if animal & workflow defined: find specific workflow over all days
%
% Output:
% recordings (struct: length = n days):
% .animal - animal
% .day - day of recording
% .recording - recording folder name(s) (time as HHMM)
% .index - index(ies) of recording within day (e.g. 3rd recording is [3])
% .workflow - Bonsai workflow(s)
% .mousecam - whether mousecam was recorded for recording(s)
% .widefield - whether widefield was recorded for recording(s)
%
% e.g. recording(4).workflow{3}: the 3rd matching workflow from the 4th day
% containing a matching recording.

%% Validate/default inputs
arguments

    animal char {mustBeNonempty}
    recording_day = {};
    workflow = {};

end

%% Initialize parameters

% Check day pattern
if ~exist('recording_day','var') || ~isempty(recording_day)
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    if ~matches(recording_day,day_pattern)
        error('Recording day (''%s'') not ''yyyy-mm-dd'' format',recording_day);
    end
end


if ~isempty(recording_day) && ~iscell(recording_day)
    recording_day = {recording_day};
end

% If recording days empty: search all days on server
if isempty(recording_day)
    % Get contents of animal path
    animal_path = fullfile(ds.locations.local_data_path,animal);
    animal_dir = dir(animal_path);

    % Find recording paths (matches day format and is folder)
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    recording_day_idx = matches({animal_dir.name},day_pattern) & [animal_dir.isdir];
    recording_day = {animal_dir(recording_day_idx).name};
end

%% Find and package matching recordings

struct_fieldnames = ...
    {'animal','day', ...
    'video_track','tetrode','ca_2p'};
recordings = cell2struct(cell(length(struct_fieldnames),0),struct_fieldnames);

if ~isempty(recording_day)

    for curr_day_idx = 1:length(recording_day)
        curr_day = recording_day{curr_day_idx};

        curr_day_path = ds.locations.filename('local',animal,curr_day);

    

        % Package info
        % (set index for day)
        recording_idx = length(recordings) + 1;

        % (basic info)
        recordings(recording_idx).animal = animal;
        recordings(recording_idx).day = curr_day;
  
        % (recording modalities - note ephys is day-, not recording-specific)
        recordings(recording_idx).video_track = ...
           any(exist(fullfile(curr_day_path,'video_track'),'dir'))
            
        recordings(recording_idx).tetrode = ...
               any(exist(fullfile(curr_day_path,'tetrode_recording'),'dir'))
        recordings(recording_idx).ca_2p = ...
            any(exist(fullfile(curr_day_path,'ca_2p'),'dir'));
    end

    %% If no recording days, error out
else
    error('No recording days: %s',animal)
end


%% If no recordings found, error out
if isempty(recordings)
    warning('No recordings found: %s [%s] [%s]',animal,strjoin(recording_day,','),strjoin(workflow,','))
end



















