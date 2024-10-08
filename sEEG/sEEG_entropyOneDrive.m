
addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/1. Recording & Stimulation of iEEG - pbelchps files/Image Reconstruction/customFcns'))
addpath(genpath('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/abelCode/ripple/Tools/neuroshare'))

datapath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/sEEG_rawData/');
patientAnatPath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/');
savePath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results');
entropyPath = ('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/individualSubFiles/');

agefile = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/PBE_agefile.csv');

%load in all the files
setfiles0 = dir([datapath, 'P*/Rest/rest*ns2*']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end


%% permutation entropy 

for j = 1:length(setfiles0)

    idvalues(j,:) = (setfiles0(j).folder(106:112));

    if ~isfile([entropyPath idvalues(j,:) '_permEntropy.csv'])

        if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
            load([savePath, idvalues(j,:), '/Rest/data.mat']);
            fiveMin = data(1:300000, :);
            dataDown = downsample(fiveMin, 2);
            montage = readtable([savePath, idvalues(j,:), '/Rest/montage', idvalues(j,:), '.xlsx']);
           
            % Identify rows where the first column has NaN in the last two rows
            rowsToRemove = isnan(montage{end-1:end, 1});

            if any(rowsToRemove)
                % Remove the identified rows
                montage(end-1:end, :) = [];
                disp('Rows removed.');
            else
                disp('No rows to remove.');
            end


            for c = 1:size(data, 2)
                [Perm(c,:), Pnorm(c,:), cPE(c,:)] = PermEn(fiveMin(:,c), 'm', 5, 'tau', 1);

            end

            permTable = array2table(Perm);
            pnormTable = array2table(Pnorm);
            cPETable = array2table(cPE);

            subjectTable = horzcat(montage, permTable, pnormTable, cPETable);
            
            
            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];

          
            subjectSavePath = [entropyPath idvalues(j,:) '_permEntropy.csv'];
            writetable(subjectTable, subjectSavePath)

        end

    end

    clear Perm
    clear permTable
    clear pnormTable
    clear cPETable
    clear Pnorm
    clear cPE

end


%% Sample entropy 

for j = 1:length(setfiles0)

    idvalues(j,:) = (setfiles0(j).folder(106:112));

    if ~isfile([entropyPath idvalues(j,:) '_permEntropy.csv'])

        if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
            load([savePath, idvalues(j,:), '/Rest/data.mat']);
            fiveMin = data(1:300000, :);
            oneMin = data(1:60000,:);
            dataDown = downsample(fiveMin, 2);
            montage = readtable([savePath, idvalues(j,:), '/Rest/montage', idvalues(j,:), '.xlsx']);
           
            % Identify rows where the first column has NaN in the last two rows
            rowsToRemove = isnan(montage{end-1:end, 1});

            if any(rowsToRemove)
                % Remove the identified rows
                montage(end-1:end, :) = [];
                disp('Rows removed.');
            else
                disp('No rows to remove.');
            end


            for c = 1:size(data, 2)
                [Samp(c,:)] = SampEn(oneMin(:,c), 'm', 5, 'tau', 1);

            end

            SampTable = array2table(Samp);
         
            subjectTable = horzcat(montage, SampTable);
            
            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];

          
            subjectSavePath = [entropyPath idvalues(j,:) '_sampEntropy.csv'];
            writetable(subjectTable, subjectSavePath)

        end

    end

    clear Perm
    clear permTable
    clear pnormTable
    clear cPETable
    clear Pnorm
    clear cPE

end




%% multiscale entropy 
parpool('local', 20);

for j = 1:length(setfiles0)

    idvalues(j,:) = (setfiles0(j).folder(80:86));
    pt = (setfiles0(j).folder(80:86));
    load([datapath, idvalues(j,:), '/Rest/data.mat']);

    chans = info.chanNames;
    % rois = chan2roi(pt, chans);
    % rois(:,6) = chans;

    if ~isfile([entropyPath idvalues(j,:) '_5mins_MultiScaleEntropy20_referenced.csv'])

        if isfile([datapath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
           
            load([datapath, idvalues(j,:), '/Rest/data.mat']);
            dataCar = data - mean(data,2);
            % oneMin = dataCar(1:60000, :);
            fiveMin = dataCar(1:300000, :);
            dataDown = downsample(fiveMin,7);
            montage = readtable([datapath, idvalues(j,:), '/montage', idvalues(j,:), '.xlsx']);
           
            % Identify rows where the first column has NaN in the last two rows
            rowsToRemove = isnan(montage{end-1:end, 1});

            if any(rowsToRemove)
                % Remove the identified rows
                montage(end-1:end, :) = [];
                disp('Rows removed.');
            else
                disp('No rows to remove.');
            end


            parfor c = 1:size(data, 2)

                Mobj = MSobject("SampEn");
                [MSx(c,:), Ci(:,c)] = MSEn(dataDown(:,c), Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);

            end

            MSxTable = array2table(MSx);
            CiTable = array2table(Ci');

            subjectTable = horzcat(montage, MSxTable, CiTable);
            
            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];

          
            subjectSavePath = [entropyPath idvalues(j,:) '_5mins_MultiScaleEntropy20_referenced.csv'];
            writetable(subjectTable, subjectSavePath)

        end

    end

    clear Mobj
    clear MSx
    clear Ci

end
 delete(gcp);


%% Spectral Entropy 


for j = 1:length(setfiles0)

    idvalues(j,:) = (setfiles0(j).folder(106:112));

    if ~isfile([entropyPath idvalues(j,:) '_spectralEntropy_referenced.csv'])

        if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
            load([savePath, idvalues(j,:), '/Rest/data.mat']);
            dataCar = data - mean(data,2);
            oneMin = dataCar(1:60000, :);
            % fiveMin = data(1:300000, :);
            dataDown = downsample(oneMin, 2);
            montage = readtable([savePath, idvalues(j,:), '/Rest/montage', idvalues(j,:), '.xlsx']);
           
            % Identify rows where the first column has NaN in the last two rows
            rowsToRemove = isnan(montage{end-1:end, 1});

            if any(rowsToRemove)
                % Remove the identified rows
                montage(end-1:end, :) = [];
                disp('Rows removed.');
            else
                disp('No rows to remove.');
            end

    
            for c = 1:size(data, 2)

               [Spec(c,:), highGammaBandEn(c,:)] = SpecEn(dataDown(:,c), 'N', 512, 'Freqs', [.312, 1], 'Logx', exp(1), 'Norm' , true);
               [Spec(c,:), gammaBandEn(c,:)] = SpecEn(dataDown(:,c), 'N', 512, 'Freqs', [.12, .312], 'Logx', exp(1), 'Norm' , true);
               [Spec(c,:), betaBandEn(c,:)] = SpecEn(dataDown(:,c), 'N', 512, 'Freqs', [.05, .1], 'Logx', exp(1), 'Norm' , true);
               [Spec(c,:), alphaBandEn(c,:)] = SpecEn(dataDown(:,c), 'N', 512, 'Freqs', [.03, .05], 'Logx', exp(1), 'Norm' , true);
               [Spec(c,:), thetaBandEn(c,:)] = SpecEn(dataDown(:,c), 'N', 512, 'Freqs', [.01, .03], 'Logx', exp(1), 'Norm' , true);

            end



            highGammaTable = array2table(highGammaBandEn);
            gammaBandEnTable = array2table(gammaBandEn);
            betaBandEnTable = array2table(betaBandEn);
            alphaBandEnTable = array2table(alphaBandEn);
            thetaBandEnTable = array2table(thetaBandEn);
            specTable = array2table(Spec);


            subjectTable = horzcat(montage, highGammaTable, gammaBandEnTable, betaBandEnTable, alphaBandEnTable, thetaBandEnTable, specTable);
            
            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];

          
            subjectSavePath = [entropyPath idvalues(j,:) '_spectralEntropy_referenced.csv'];
            writetable(subjectTable, subjectSavePath)

        end

    end

    clear highGammaBandEn
    clear gammaBandEn
    clear betaBandEn
    clear alphaBandEn
    clear thetaBandEn
    clear subjectTable
    clear Spec


end