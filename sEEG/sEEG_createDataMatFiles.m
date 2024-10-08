
addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/seeg/Image Reconstruction/customFcns'))
addpath(genpath('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/abelCode/ripple/Tools/neuroshare'))

datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/sEEG_rawData/');
patientAnatPath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive - University of Pittsburgh/seeg/Image Reconstruction/patients/');
savePath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/sEEG_rawData/');

%load in all the files
setfiles0 = dir([datapath, 'P*/Rest/rest*ns2*']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end

for j = 1:length(setfiles0)
    idvalues(j,:) = (setfiles0(j).folder(80:86));
    disp(idvalues(j,:))
    if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
        continue
    else % no .mat file for patient. we need to make it
        [nsResult,hFile] = ns_OpenFile(fullfile([setfiles0(j).folder,'/', setfiles0(j).name]));
        mtgTbl = readtable(fullfile([datapath, idvalues(j,:), '/montage', idvalues(j,:),'.xlsx']));
        [chans,isort] = sort(mtgTbl{~isnan(mtgTbl{:,1}),2});
        info.chanNames = chans;
        % get indices of all channels labeled 'lfp #'
        labels = {hFile.Entity.Label};
        iLfp = find(cellfun(@(x) contains(x,'lfp '),labels));
        % extract block of all data
        [nsResult,dataBlock] = ns_GetAnalogDataBlock(hFile,iLfp,1,hFile.Entity(iLfp(1)).Count);
        data = dataBlock(:,isort);
        
        save(fullfile([savePath, idvalues(j,:), '/Rest/data.mat']),'info','data','-v7.3')
        
        
    end
    
end

for j = 1:length(setfiles0)
    idvalues(j,:) = (setfiles0(j).folder(99:105));
    disp(idvalues(j,:))
    if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient, if there is lets make an roi file
        
            anatFolder = ([patientAnatPath, idvalues(j,:)]);
            
            if isfile([patientAnatPath, idvalues(j,:), '/anat.mat']) % is there a anat file
                load([anatFolder, '/anat.mat']); % load anat file
                load([savePath, idvalues(j,:), '/data.mat']); % load data file
                
                roi = chan2roi(idvalues(j,:),info.chanNames); % run chan2roi function
                
                save([savePath, idvalues(j,:), '/roi.mat'], 'roi')
                writecell(roi,[savePath, idvalues(j,:), '/roi.csv']);
            else %if there is no anat file, we have to break the loop, cause we cant make an roi file
                continue
                
            end
      
    else % if there isnt a data.mat file, break the loop and move on
        continue
    end

end



