
%% set initial paths
addpath(genpath(('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1')));
addpath(genpath(('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/code')));
datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/Resting_State/AfterWhole/ICAwholeClean_homogenize');
outpath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/');

%% load in all the data files
setfiles0 = dir([datapath,'/*icapru*.set']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1 : length(setfiles0)
    idvalues{j} = (setfiles0(j).name(1:14));
end

numSubj = length(idvalues);
closedeyes_codes = {'16129', '15261', '0', '15361'};

i = 1;
for i = 1:numSubj
    disp(i);
    subject = idvalues{i};
    inputfile = setfiles{i};

    savePath = [outpath subject '_acw_eyesClosed_newACWscript.csv'];

    if ~exist(savePath, 'file')
        fprintf('The file %s does not exist, running ACW code\n', savePath)
        % open eeglab
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        try
            EEG = pop_loadset(inputfile); % load in eeg file
        catch
            disp(['subject missing fdt: ' subject])
            continue;
        end
    end

    try
        %EEGopeneyes = pop_rmdat(EEG, {'16130', '15362','1', '15362'},[0 4] ,0);
        
       EEGclosedeyes = pop_rmdat(EEG, {'16129', '15261','0', '15361'},[0 4] ,0);
       
       % Check if the remaining events match the desired event codes
        remaining_events = {EEGclosedeyes.event.type}; % Get the event types in EEGclosedeyes

        % Find the intersection between remaining events and the desired event codes
        matching_events = intersect(remaining_events, closedeyes_codes);

        % If no matching events are found, skip to the next subject
        if isempty(matching_events)
            disp(['No matching events found for subject ' num2str(subject) ', skipping...']);
            continue; % Skip to the next subject
        end
    
        for c = 1:size(EEG.data,1)
            ACWout{c,1} = EEGclosedeyes.chanlocs(c).labels;            
            [ACWout{c,2}] = ACW_estimation(EEGclosedeyes.data(c,:),EEGclosedeyes.srate,20);
            [ACWout{c,3}, ACWout{c,4}, acf, lags] = acw(EEGclosedeyes.data(c,:),EEGclosedeyes.srate, 0);
        end

        % Create a new column with subject ID repeated for every row
        subjectIDColumn = repmat(idvalues(:,i), size(ACWout, 1),1);

        % Add the new column to the existing table
        subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), ACWout];
        subjectTable.Properties.VariableNames{2} = 'Channel';
        subjectTable.Properties.VariableNames{3} = 'ACW_old';
        subjectTable.Properties.VariableNames{4} = 'ACW_0';
        subjectTable.Properties.VariableNames{5} = 'ACW_50';

        subjectSavePath = [outpath subject '_acw_eyesClosed_newACWscript.csv'];
        writetable(subjectTable, subjectSavePath)
        
    catch
        disp(['subject not running through: ' subject])
        continue;

    end

end
