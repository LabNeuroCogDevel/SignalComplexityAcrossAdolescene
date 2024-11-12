
function [subjectTable] = run_Entropy_PSC(inputfile, task, epoch, lengthValue, varargin)

addpath(genpath(('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/resources/eeglab2022.1')));
addpath(('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/resources/fieldtrip-20220104'))
ft_defaults

if task == 'MGS'
    resultFolder = 'MGS_Entropy';
elseif task == 'rest'
    resultFolder = 'entropy';
end

resultPath = ('/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/Results');

savePath = [resultPath '/' resultFolder '/individual_subject_files/'];

EEG = pop_loadset(inputfile); % load in eeg file
[d, currentName, ext ] = fileparts(inputfile);


if strcmp(task, 'MGS') && strcmp(epoch, 'delay')
    if ~isfile([savePath idvalues(i,:) '_MultiScaleEntropy_delay' num2str(lengthValue) '.csv'])

        if size(EEG.data,1) > 64
            EEG = pop_select( EEG,'nochannel',{'EX3' 'EX4' 'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});
        end

        if ~isfield(EEG.event(1), 'seconds')
            for e = 1:length(EEG.event)
                %find event latencies in seconds
                EEG.event(e).seconds = ([EEG.event(e).latency]-1)/EEG.srate;
            end
        end

        for e = 1:length(EEG.event)
            %find event latencies in seconds
            if e == 1
                EEG.event(e).duration = EEG.event(e).seconds;
            elseif e == length(EEG.event)
                EEG.event(e).duration = 2;
            else
                EEG.event(e).duration = round(EEG.event(e+1).seconds - EEG.event(e).seconds);
            end
        end

        %select delay epochs based on length
        for e = 1:length(EEG.event)
            if num2str(EEG.event(e).type) == '4'
                EEG.event(e).type = ['4_' num2str((EEG.event(e).duration))];
            end
        end

        try

            if lengthValue == 6
                delayEEG = pop_selectevent( EEG, 'type',{'4_6'},'deleteevents','on');
                % Find indices of boundary events
                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];

                delayEEG = pop_rmdat(delayEEG, {'4_6'},[0 5.98] ,0);

                % Find indices of boundary events
                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];

            elseif lengthValue == 8
                delayEEG = pop_selectevent( EEG, 'type',{'4_8'},'deleteevents','on');
                % Find indices of boundary events
                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];

                delayEEG = pop_rmdat(delayEEG, {'4_8'},[0 7.98] ,0);

                % Find indices of boundary events
                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];


            elseif lengthValue == 10
                delayEEG = pop_selectevent( EEG, 'type',{'4_10'},'deleteevents','on');
                % Find indices of boundary events
                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];

                delayEEG = pop_rmdat(delayEEG, {'4_10'},[0 9.98] ,0); % Find indices of boundary events

                boundary_indices = find(strcmp({delayEEG.event.type}, 'boundary'));

                % Remove boundary events from EEG.event structure
                delayEEG.event(boundary_indices) = [];

            end

            if (length(delayEEG.event) == length(EEG.event)) || (length(delayEEG.event) < 50 && lengthValue == 6) || (length(delayEEG.event) < 20 && lengthValue == 8)|| (length(delayEEG.event) < 5 && lengthValue == 10)
                return; % Move to next person in the loop if no events were removed
            end


            [subjectTable] = Calculate_EEG_Entropy_Values(delayEEG);

            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(i,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];


            subjectSavePath = [savePath idvalues(i,:) '_MultiScaleEntropy_delay' num2str(lengthValue) '.csv'];
            writetable(subjectTable, subjectSavePath);
        catch
            return;
        end


    end

elseif task == 'rest'

    if ~isfile([savePath idvalues(i,:) '_MultiScaleEntropy_eyesClosed.csv'])
        EEGclosedeyes = pop_rmdat(EEG, {'16129', '15261','0'},[0 4] ,0);

        [subjectTable] = Calculate_EEG_Entropy_Values(EEGclosedeyes);

        % Create a new column with subject ID repeated for every row
        subjectIDColumn = repmat(idvalues(i,:), size(subjectTable, 1),1);

        % Add the new column to the existing table
        subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];


        subjectSavePath = [savePath idvalues(i,:) '_MultiScaleEntropy_eyesClosed.csv'];
        writetable(subjectTable, subjectSavePath);


    end
end

