%% ***** Formal Coding V1.0 *****

% Cued responses in cued rewarded trials
% Test gap from 3, 6, 9, 12, 15 frames between movements

% The aim is to extract the behavior components, including: successful ratio, 
% time, speed, amplitude, movement ratio, etc.
% Further unsed for look at corresponding frames in imaging.
% This code works on multiple days for individual animal.

% Current code is based on CR_Image_LeverPress_MultipleDay and
% AP_leverpress_bhv.
% 
IN = 'VIP';
Animals = {'3438544-O','3438544-L','3438544-R'};

Gaps = round(3*33.33)*5; % 500 ms
Leeway = 150; % 150 ms
MovDurs = 200; % 2 s

cd(['Z:\People\Chi\WFLP_IN' filesep IN]);

for curr_gap = 1:length(Gaps)
    
    for curr_animal = 1:length(Animals)

        %% Clean
        clearvars -except Animals IN curr_animal Gaps curr_gap MovDurs MovDurs_Rewarded Leeway;
        clc;

        %% Set animal name, sample range, and threshold of correlation coefficient of movement
        % The initial format is written according to the organization of files in
        % Chi's folder on the server.
        % Ideally all changes (animal name, data, sample range) to be made are in this section.

        % *** CHANGE SECTION STARTS *** %
        % Set animal name and date
        Initial = 'CR';
        Animal = Animals{curr_animal};
        Gap = Gaps(curr_gap);
        MovDur = MovDurs(curr_gap);
        % Set sample range
        Sample_range_xsg = [-5000:40000]; % 1 ms = 10 frames in xsg;
        Sample_range_resample = [-500:4000];
        Sample_range_corr = [50:50+MovDur]; % 10 ms step
        Lever_downsamp = 100;

        % *** CHANGE SECTION ENDS *** %
        
        cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal]);
        mkdir('MovAnalysis');
        
        if ismember(Animal,{'3271320-L','3333408-O','3333408-L','3333408-R'})
            task_start_date = '171105';
        elseif ismember(Animal,{'3258531-O','3258531-L'})
            task_start_date = '171129';
        elseif ismember(Animal, {'3358884-L','3358884-R'})
            task_start_date = '171213';
        elseif ismember(Animal, {'3358883-O','3358883-L','3358883-R'})
            task_start_date = '171214';
        elseif ismember(Animal, {'3373693-O','3373693-L'})
            task_start_date = '180111';
        elseif ismember(Animal, {'3373693-R','3373693-LR'})
            task_start_date = '180117';
        elseif ismember(Animal, {'3420509-L','3420509-R'})
            task_start_date = '180206';       
        elseif ismember(Animal, {'3491479-L','3438483-L'})
            task_start_date = '180501';
        elseif ismember(Animal, {'3491479-LR'})
            task_start_date = '180503';
        elseif ismember(Animal, {'3491479-R','3438500-O','3438500-L'})
            task_start_date = '180508';
        elseif ismember(Animal, {'3438521-R','3438521-LR'})
            task_start_date = '180530';
        elseif ismember(Animal, {'3438521-O','3438521-L'})
            task_start_date = '180531';
        else
            task_start_date = [];
        end
            

        %% Load files
        % Load data path
        if ispc
            Data_path = ['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal];
        elseif ismac
            Data_path = ['/Volumes/lab/People/Chi/WFLP_IN' filesep IN filesep Initial '_' Animal];
        end
        
        %% Find out data path and dates
        % For dispatcher
        Dispatcher_path = [Data_path filesep 'Dispatcher'];
        Dispatcher_file_list = dir(Dispatcher_path); % strcuture
        Dispatcher_file_list = Dispatcher_file_list(3:end); % The first two are '.' and '..'
        % For xsg
        Xsg_path = [Data_path filesep 'LeverTrace'];
        Xsg_folder_list = dir(Xsg_path);
        Xsg_folder_list = Xsg_folder_list(3:end);
        Date_list = cellfun(@(x) x(1:6), {Xsg_folder_list.name}, 'UniformOutput' ,false);
        
        if isempty(task_start_date)
            index = 1;
        else
            index = find(ismember(Date_list,task_start_date));
            Date_list = Date_list(index:end);
        end
        
        % Load trial index
        load(['Z:\People\Chi\WFLP_IN\CheckNonconsTrials' filesep IN '_' Initial '_' Animal '_CheckNonconsTrials_Corrected.mat'], 'TrialInfo_xsg_all','SessionTerm');
        TrialInfo_xsg_all = TrialInfo_xsg_all(index:end);
        %% Load Dispatcher and xsg files for current day
        disp([Initial '_' Animal ' ' num2str(Gap)]);

        % Initialization
        CuedMov_SingleAnimal = cell(size(Date_list));

        % Basically go throuh trial by trial
        for curr_day = 1:length(Date_list)
            
            licking_threshold = 4;
            Date = Date_list{curr_day};
            % Should only have one filename
            Curr_Dispatcher_filname = {Dispatcher_file_list(cellfun(@(x) ~isempty(strfind(x, Date)), {Dispatcher_file_list.name})).name};
            if length(Curr_Dispatcher_filname) > 1
                warning(['Incorrect Dispatcher file saved for ' Animal ' on ' Date '. Please check.'])
                return
            end
            if isempty(Curr_Dispatcher_filname)
                warning(['No Dispatcher file saved for ' Animal ' on ' Date '.']);
                continue
            end
            warning off
            load([Dispatcher_path filesep Curr_Dispatcher_filname{1}]); % Because cd is not in the Dispatcher folder.
            warning on

            % Get general trial information that will be used for all xsg files
            % Extract trial information from dispatcher
            TrialInfo_Dispatcher = saved_history.ProtocolsSection_parsed_events;
            CatchInfo_Dispatcher = cell2mat(saved_history.CatchSection_catch_trial);
            % Total trial num
                
%             if strcmp(IN,'ChAT_Cas')
%                 TrialNum_Total = SessionTerm(curr_day); 
%             elseif strcmp(IN,'VIP_ArchT') && strcmp(Animal,'4412540-R')
%                 TrialNum_Total = SessionTerm(curr_day); 
            if exist('SessionTerm')
                TrialNum_Total = SessionTerm(curr_day); 
            else
                TrialNum_Total = length(TrialInfo_Dispatcher);
            end
            TrialInfo_Dispatcher = TrialInfo_Dispatcher(1:TrialNum_Total);
            % Get rewarded trials
            Rewarded_trial = find(cellfun(@(x) ~isempty(x.states.reward), TrialInfo_Dispatcher));

            %% Load xsg file list
            Xsg_curr_path = [Data_path filesep 'LeverTrace' filesep Date];
            Xsg_file_list = dir(Xsg_curr_path);
            % Should have multiple filenames, use cell, otherwise will only save the
            % first one.
            Xsg_filename_list = {Xsg_file_list(cellfun(@(x) ~isempty(strfind(x, 'CR0001AAAA')), {Xsg_file_list.name})).name};
            if isempty(Xsg_filename_list)
                warning(['No xsg file saved for ' Animal ' on ' Date '.']);
                CuedMov.Cued_MovOnset_Info = {};
                CuedMov.LeverTrace = {};
                CuedMov.Cued_MovOnset_Info_All = [];
                CuedMov.LeverTrace_All = [];
                CuedMov.LeverTrace_First = [];
                CuedMov.LeverTrace_Reward = [];
                CuedMov.LeverTrace_All_downsample = [];
                CuedMov.LeverTrace_First_downsample = [];
                CuedMov.LeverTrace_Reward_downsample = [];
                CuedMov.Median_Movduration = nan;
                CuedMov.Median_Movduration_First = nan;
                CuedMov.Median_Movduration_Reward = nan;
                CuedMov.CR = length(Rewarded_trial)/TrialNum_Total;
                CuedMov_SingleAnimal{curr_day} = CuedMov;
                continue
            end
            Xsg_filename_list = sort(Xsg_filename_list); % from 1st to end

            % Initialization for each day
            TrialInfo_xsg = TrialInfo_xsg_all{curr_day};
            CuedMov.Cued_MovOnset_Info = cell(size(Xsg_filename_list));
            CuedMov.LeverTrace = cell(size(Xsg_filename_list));
            CuedMov.Cued_MovOnset_Info_All = [];
            CuedMov.LeverTrace_All = [];
            CuedMov.LeverTrace_First = [];
            CuedMov.LeverTrace_Reward = [];
            CuedMov.LeverTrace_All_downsample = [];
            CuedMov.LeverTrace_First_downsample = [];
            CuedMov.LeverTrace_Reward_downsample = [];
            CuedMov.Median_Movduration = nan;
            CuedMov.Median_Movduration_First = nan;
            CuedMov.Median_Movduration_Reward = nan;
            CuedMov.CR = length(Rewarded_trial)/TrialNum_Total;
            StiTime.CueOn{curr_day} = cell(size(Xsg_filename_list));
            StiTime.CueOff{curr_day} = cell(size(Xsg_filename_list));
            StiTime.RewardOn{curr_day} = cell(size(Xsg_filename_list));
            StiTime.RewardOff{curr_day} = cell(size(Xsg_filename_list));
            StiTime.PunishOn{curr_day} = cell(size(Xsg_filename_list));
            StiTime.PunishOff{curr_day} = cell(size(Xsg_filename_list));
            %% Load xsg file sequentially

            % User interaction
            disp(['Running ' Animal ' on ' Date]);

            for curr_xsg = 1:length(Xsg_filename_list)
                if strcmp(Animal,'3233232-R') && strcmp(Date,'170723') && curr_xsg == 1
                    continue
                end
                % Ignore if no trials exist
                if isempty(TrialInfo_xsg{curr_xsg})
                    disp(['No Bitecode in ' Xsg_filename_list{curr_xsg} '.'])
                    continue
                end
                % Loading
                load([Xsg_curr_path filesep Xsg_filename_list{curr_xsg}],'-MAT'); % Because cd is not in the LeverTrace folder.
                % Get xsg sample rate
                Sample_rate_xsg = header.acquirer.acquirer.sampleRate; % Should be 10000
                % Get recording of bitcode, lever trace and licking
                Lever_trace_xsg = data.acquirer.trace_2;
                Licking_xsg = data.acquirer.trace_3;
                if isfield(data.acquirer, 'trace_4')
                    PupilFrameStamp_xsg = data.acquirer.trace_4;
                end
               
                % For bad linux machine: trial index error reading
                consec_trials = diff([0; TrialInfo_xsg{curr_xsg}(:,2)]);
                if sum(consec_trials(2:end)~=1) > 1
                    disp('Non-consecutive trials');
                end
                
                % Ignore if no trials exist after correcting
                if isempty(TrialInfo_xsg{curr_xsg})
                    disp(['No Bitecode in ' Xsg_filename_list{curr_xsg} ' after correction.'])
                    continue
                end
                % For bad linux machine: Sometimes the last trial index of last xsg will exceeds dispatcher if
                % we stop between bitcode and cue, very rare I guess
                if TrialInfo_xsg{curr_xsg}(end,2) > TrialNum_Total
                    TrialInfo_xsg{curr_xsg}(end,:) = [];
                end
                % Ignore if no trials exist after correcting
                if isempty(TrialInfo_xsg{curr_xsg})
                    disp(['No Bitecode in ' Xsg_filename_list{curr_xsg} ' after correction.'])
                    continue
                end

                % Get movment start and stop from lever trace
                [lever_active lever_force_resample lever_force_smooth ...
                    lever_velocity_envelope_smooth] ...
                    = CR_parseLeverMovement_Updated(Lever_trace_xsg, Gap, Leeway); % 1K Hz
                Lever_active_switch = diff([0;lever_active]); % 1K Hz                
                Lever_Active{curr_day}{curr_xsg,1} = lever_active;
                Lever_Force_Resample{curr_day}{curr_xsg,1} = lever_force_resample;
                Lever_MovOnset{curr_day}{curr_xsg,1} = find(Lever_active_switch == 1);
                
                [Licking_active, Licking_bout, Licking_bout_switch] = CR_GetLickingState(Licking_xsg,licking_threshold);
                Licking_Active{curr_day}{curr_xsg,1} = Licking_active;
                Licking_Bout{curr_day}{curr_xsg,1} = Licking_bout;
                Licking_Onset{curr_day}{curr_xsg,1} = find(Licking_bout_switch == 1);
                
                % Get frame timestamp for pupil imaging
                if isfield(data.acquirer, 'trace_4')
                    PupilFrameStamp_xsg = PupilFrameStamp_xsg > 4;
                    PupilFrameStamp_xsg = diff([0;PupilFrameStamp_xsg]);
                    PupilFrameStamp{curr_day}{curr_xsg,1} = round(find(PupilFrameStamp_xsg == 1)./10)+1;
                else
                    PupilFrameStamp{curr_day}{curr_xsg,1} = [];
                end

                MovOnset_curr_xsg = [];
                LeverTrace_curr_xsg = [];

                for curr_trial = 1:size(TrialInfo_xsg{curr_xsg},1)
                    
                    curr_trial_index = TrialInfo_xsg{curr_xsg}(curr_trial,2);
                    % Ignore if unrewarded trials
                    if curr_trial_index > TrialNum_Total
                        continue
                    end
                    
                    % Get stimulus time
                    Bitcode_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.bitcode(1);
                    StiTime.CueOn{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.cue(1)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                    StiTime.CueOff{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.cue(2)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                    if ~isempty(TrialInfo_Dispatcher{curr_trial_index,1}.states.reward)
                        StiTime.RewardOn{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.reward(1)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                        StiTime.RewardOff{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.reward(2)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                    else
                        StiTime.RewardOn{curr_day}{curr_xsg}(curr_trial,1) = nan;
                        StiTime.RewardOff{curr_day}{curr_xsg}(curr_trial,1) = nan;
                    end
                    if ~isempty(TrialInfo_Dispatcher{curr_trial_index,1}.states.end_cue_punish)
                        StiTime.PunishOn{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.end_cue_punish(1)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                        StiTime.PunishOff{curr_day}{curr_xsg}(curr_trial,1) = TrialInfo_Dispatcher{curr_trial_index,1}.states.end_cue_punish(2)-Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                    else
                        StiTime.PunishOn{curr_day}{curr_xsg}(curr_trial,1) = nan;
                        StiTime.PunishOff{curr_day}{curr_xsg}(curr_trial,1) = nan;
                    end
                                                                
                    
                    if ismember(curr_trial_index, Rewarded_trial)
                        Bitcode_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.bitcode(1);
                        % Cue time
                        Cue_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.cue(1);
                        Cue_xsg = Cue_dispatcher - Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                        % Reward time
                        Reward_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.reward(1);
                        Reward_xsg = Reward_dispatcher - Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                        % Refine, within 300 sec
                        if Cue_xsg>300 || Reward_xsg>300
                            continue
                        end
                        % Find out it's cued or not
                        if ~isempty(find(Lever_active_switch(1:round(Reward_xsg*1000))==1,1,'last'))
                            Rewarded_MovOnset = find(Lever_active_switch(1:round(Reward_xsg*1000))==1,1,'last')/1000; % sec
                            if Rewarded_MovOnset >= Cue_xsg && sum(lever_active((round(Cue_xsg*1000)-100):round(Cue_xsg*1000)))==0 && Rewarded_MovOnset+Sample_range_resample(end)/1000<=300 % Collect trace
                                Temp_MovOnset_info = [];
                                Temp_MovTrace_curr_trial = [];
                                MovOnset_curr_trial = [];
                                MovOnset_curr_trial = find(Lever_active_switch(round(Cue_xsg*1000):round(Reward_xsg*1000))==1) + round(Cue_xsg*1000) - 1;
                                for curr_movonset = 1:length(MovOnset_curr_trial)
                                    Temp_MovOnset_info(curr_movonset,1) = curr_trial_index;
                                    Temp_MovOnset_info(curr_movonset,2) = MovOnset_curr_trial(curr_movonset)./1000; % sec
                                    % Get movement duration
                                    Temp_MovOnset_info(curr_movonset,3) = find(Lever_active_switch(MovOnset_curr_trial(curr_movonset):end)==-1,1,'first')/1000;
                                    if isnan(Temp_MovOnset_info(curr_movonset,3))
                                        Temp_MovOnset_info(curr_movonset,3) = 300 - Temp_MovOnset_info(curr_movonset,2);
                                    end
                                    Temp_MovOnset_info(curr_movonset,4) = curr_movonset;
                                    if MovOnset_curr_trial(curr_movonset)/1000 == Rewarded_MovOnset
                                        Temp_MovOnset_info(curr_movonset,5) = true;
                                    else
                                        Temp_MovOnset_info(curr_movonset,5) = false;
                                    end
                                    Temp_MovOnset_info(curr_movonset,6) = Cue_xsg;
                                    Temp_MovOnset_info(curr_movonset,7) = Reward_xsg;
                                    Temp_MovOnset_info(curr_movonset,8) = CatchInfo_Dispatcher(curr_trial_index);
                                    % Licking onset
                                    if ~isempty(find(Licking_Onset{curr_day}{curr_xsg,1}./1000>Reward_xsg,1,'first'))
                                        Temp_MovOnset_info(curr_movonset,9) = Licking_Onset{curr_day}{curr_xsg,1}(find(Licking_Onset{curr_day}{curr_xsg,1}./1000>Reward_xsg,1,'first'))./1000;
                                    else
                                        Temp_MovOnset_info(curr_movonset,9) = nan;
                                    end
                                    Temp_MovTrace_curr_trial = horzcat(Temp_MovTrace_curr_trial,Lever_trace_xsg(round(MovOnset_curr_trial(curr_movonset)*10+Sample_range_xsg)));
                                end
                                MovOnset_curr_xsg = vertcat(MovOnset_curr_xsg,Temp_MovOnset_info);
                                LeverTrace_curr_xsg = horzcat(LeverTrace_curr_xsg,Temp_MovTrace_curr_trial);
                            end
                        end
                    end
                end
                disp(['Finish ' Xsg_filename_list{curr_xsg}]);
                CuedMov.Cued_MovOnset_Info{curr_xsg} = MovOnset_curr_xsg;
                CuedMov.LeverTrace{curr_xsg} = LeverTrace_curr_xsg;
            end
            CuedMov.Cued_MovOnset_Info_All = cell2mat(CuedMov.Cued_MovOnset_Info');
            CuedMov.LeverTrace_All = cell2mat(CuedMov.LeverTrace);
            if ~isempty(CuedMov.LeverTrace_All)
                Temp_trace = resample(CuedMov.LeverTrace_All,1,Lever_downsamp);
                CuedMov.LeverTrace_All_downsample = Temp_trace(Sample_range_corr,:);
            else
                CuedMov.LeverTrace_All_downsample = [];
            end
            % First movement
            if ~isempty(CuedMov.LeverTrace_All)
                First_index = find(CuedMov.Cued_MovOnset_Info_All(:,4)==1);
                CuedMov.LeverTrace_First = CuedMov.LeverTrace_All(:,First_index);
                Temp_trace = resample(CuedMov.LeverTrace_First,1,Lever_downsamp);
                CuedMov.LeverTrace_First_downsample = Temp_trace(Sample_range_corr,:);
            else
                CuedMov.LeverTrace_First = [];
                CuedMov.LeverTrace_First_downsample = [];
            end
            % Reward movemen
            if ~isempty(CuedMov.LeverTrace_All)
                Reward_index = logical(CuedMov.Cued_MovOnset_Info_All(:,5));
                CuedMov.LeverTrace_Reward = CuedMov.LeverTrace_All(:,Reward_index);
                Temp_trace = resample(CuedMov.LeverTrace_Reward,1,Lever_downsamp);
                CuedMov.LeverTrace_Reward_downsample = Temp_trace(Sample_range_corr,:);
            else
                CuedMov.LeverTrace_Reward = [];
                CuedMov.LeverTrace_Reward_downsample = [];
            end

            if ~isempty(CuedMov.Cued_MovOnset_Info_All)
                CuedMov.Median_Movduration = median(CuedMov.Cued_MovOnset_Info_All(:,3));
                CuedMov.Median_Movduration_Reward = median(CuedMov.Cued_MovOnset_Info_All(Reward_index,3));       
                CuedMov.Median_Movduration_First = median(CuedMov.Cued_MovOnset_Info_All(First_index,3));
            end
            CuedMov_SingleAnimal{curr_day} = CuedMov;
            disp('Finish all xsg files');
            disp(['Finish ' Date]);
        end

        %% Calculate Correlation Matrix
        disp('Calculating correlation...')
        % Trail_Trial_Correlation Matrix
        % All traces
        Trial_Trial_Corr_All = nan(length(CuedMov_SingleAnimal));
        for curr_day_1 = 1:length(CuedMov_SingleAnimal)
            if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample)
                for curr_day_2 = 1:length(CuedMov_SingleAnimal)
                    if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_All_downsample) 
                        if curr_day_1 == curr_day_2
                            Trial_Trial_Corr_All(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample);
                        else
                            temp_matrix = [];
                            temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample, ...
                                CuedMov_SingleAnimal{curr_day_2}.LeverTrace_All_downsample);
                            temp_corr_matrix = corrcoef(temp_matrix);
                            trialnum_1 = size(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample,2);
                            trialnum_2 = size(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_All_downsample,2);
                            temp_corr_matrix = temp_corr_matrix(1:trialnum_1,trialnum_1+1:end);                        
                            Trial_Trial_Corr_All(curr_day_1,curr_day_2) = nanmedian(temp_corr_matrix(:));
                        end
                    end
                end
            end
        end

        % First traces
        Trial_Trial_Corr_First = nan(length(CuedMov_SingleAnimal));
        for curr_day_1 = 1:length(CuedMov_SingleAnimal)
            if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample)
                for curr_day_2  = 1:length(CuedMov_SingleAnimal)
                    if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_First_downsample) 
                        if curr_day_1 == curr_day_2
                            Trial_Trial_Corr_First(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample);
                        else
                            temp_matrix = [];
                            temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample, ...
                                CuedMov_SingleAnimal{curr_day_2}.LeverTrace_First_downsample);
                            temp_corr_matrix = corrcoef(temp_matrix);
                            trialnum_1 = size(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample,2);
                            trialnum_2 = size(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_First_downsample,2);
                            temp_corr_matrix = temp_corr_matrix(1:trialnum_1,trialnum_1+1:end);                        
                            Trial_Trial_Corr_First(curr_day_1,curr_day_2) = nanmedian(temp_corr_matrix(:));
                        end
                    end
                end
            end
        end

        % Reward traces
        Trial_Trial_Corr_Reward = nan(length(CuedMov_SingleAnimal));
        for curr_day_1 = 1:length(CuedMov_SingleAnimal)
            if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample)
                for curr_day_2  = 1:length(CuedMov_SingleAnimal)
                    if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_Reward_downsample) 
                        if curr_day_1 == curr_day_2
                            Trial_Trial_Corr_Reward(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample);
                        else
                            temp_matrix = [];
                            temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample, ...
                                CuedMov_SingleAnimal{curr_day_2}.LeverTrace_Reward_downsample);
                            temp_corr_matrix = corrcoef(temp_matrix);
                            trialnum_1 = size(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample,2);
                            trialnum_2 = size(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_Reward_downsample,2);
                            temp_corr_matrix = temp_corr_matrix(1:trialnum_1,trialnum_1+1:end);                        
                            Trial_Trial_Corr_Reward(curr_day_1,curr_day_2) = nanmedian(temp_corr_matrix(:));
                        end
                    end
                end
            end
        end
        disp('Saving...');


        %% Save

        % User interaction
        disp('Saving...');

        % Save
        cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'MovAnalysis']);

        curr_filename_lever = [Initial '_' Animal '_CuedMov'];
        save(curr_filename_lever, 'CuedMov_SingleAnimal','Trial_Trial_Corr_All','Trial_Trial_Corr_First','Trial_Trial_Corr_Reward','-v7.3');
        curr_filename_Traces = [Initial '_' Animal '_Traces.mat'];
        if exist(curr_filename_Traces,'file')
            save(curr_filename_Traces, 'Lever_Active','Lever_Force_Resample','Lever_MovOnset','PupilFrameStamp','Licking_Active','Licking_Bout','Licking_Onset','StiTime','-append');
        else
            save(curr_filename_Traces, 'Lever_Active','Lever_Force_Resample','Lever_MovOnset','PupilFrameStamp','Licking_Active','Licking_Bout','Licking_Onset','StiTime','-v7.3');
        end
        % User interaction
        disp('Saving done');
        disp(['Finish Gap ' num2str(Gap) 'ms.']);
    end
    disp('Finish All animals');
    
    CONTINUE = false;
    if CONTINUE
        %% Caculate mean rewarded movement duration
        All.MovDuration = cell(1,9);

        for curr_animal = 1:length(Animals)
            Animal = Animals{curr_animal};

            disp(Animal);

            load(['C:\Lab\Projects\WideFieldLeverPress\GapTest' filesep 'Gap' num2str(Gap) 'ms' filesep Initial '_' Animal '_CuedMov'])
            days = length(CuedMov_SingleAnimal);

            for curr_session = 1:days
                if ~isempty(CuedMov_SingleAnimal{curr_session})
                    All.MovDuration{:,curr_animal} = vertcat(All.MovDuration{:,curr_animal},CuedMov_SingleAnimal{curr_session}.Cued_MovOnset_Info_All(:,3));
                end
            end
        end

        All.Mean_MovDuration_Single = cellfun(@nanmean,All.MovDuration);
        All.Mean_MovDuration = nanmean(cell2mat(All.MovDuration'));
        All.SEM_MovDuration = nanstd(cell2mat(All.MovDuration'))/sqrt(length(cell2mat(All.MovDuration')));

        figure(10);
        hold on
        plot(ones(1,9).*curr_gap,All.Mean_MovDuration_Single,'.','color',[0.5 0.5 0.5]);
        hold on
        plot(curr_gap,All.Mean_MovDuration,'*r');
        errorbar(curr_gap,All.Mean_MovDuration,All.SEM_MovDuration,'r'); 
    end
    
end
disp('Finish all gap values.');

