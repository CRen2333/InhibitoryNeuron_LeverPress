%% ***** Formal Coding V1.0 *****

% Purely for frame index extracting

% For IN animals, personaly I want to explore:
% 1) Cue evoked response (may not see since AC is out of current: FOV), ignore cue in the middle of movement:
% including:
%   i) responded vs non-responded
%   ii) in responded trials: rewarded or not (I guess this doesn't matter)
% 2) Rewarded movement (aligned to movement onset, calculate the delay between crossing the lower threshold and reward time, also licking on/off):
%   i) cued rewarded vs non-cued rewarded vs spontaneous initiated qualifeid movement
%   in ITI
% 3) Reward evoked response (aligned to reward)
%   i) rewarded vs spontaneous initiated qualifeid movement

% Data format:
% Column: | 1           | 2         | 3     | 4          | 5             |
% Data:   | Trial Index | Cued time | Cued? | Responded? | Reaction time |
% Column: | 6           | 7                 | 8                 |
% Data:   | Reward time | RewardedMov onset | Movement Duration |
% Column: | 9                         |
% Data:   | Crossing Higher threshold |
% Column: | 10                       | 11         | 12          | 13
% Data:   | Crossing Lower threshold | Lick onset | Lick offset | Catch

% ITI Data format: Initiated and terminated in ITI
% Column: | 1           | 2         | 3                 | 4                         |
% Data:   | Trial Index | Mov onset | Movement Duration | Crossing Higher threshold |
% Column: | 5                        | 6          | 7           | 
% Data:   | Crossing Lower threshold | Lick onset | Lick offset |

%%
IN = 'VIP';
Animals = {'3438544-O','3438544-L','3438544-R'};
Imaging_Day = [1,3,5,7,9,11,13,15,17,19,21];

Gap = round(3*33.33)*5;
Leeway = 150;

cd(['Z:\People\Chi\WFLP_IN' filesep IN]);
    
for curr_animal = 1:length(Animals)

    %% Clean
    clearvars -except Initial Animals IN curr_animal Gap Leeway Imaging_Day;
    clc;
    
    % Animal
    Animal = Animals{curr_animal};
        
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal]);
    mkdir('MovAnalysis');
    
    % Thresholds and # of imaging blocks
    switch Animal
        case '3271320-L'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [5,2,2,2,2,2,2,2,3,2,2];
            task_start_date = '171105';
        case '3333408-L'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,2,3,2,2,2,2,2,3,2,2];
            task_start_date = '171105';
        case '3333408-O'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,3,2,3,2,2,2,2,2,2,2];
            task_start_date = '171105';
        case '3333408-R'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,3,2,3,2,2,2,2,4,2,2];
            task_start_date = '171105';
        case '3258531-O'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.65,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,4,2,3,3,3,2,2,3,2,2,4];
            task_start_date = '171129';
        case '3258531-L'
            Threshold_Low = [-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.65,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [5,5,3,3,3,2,2,2,2,2,2,4];
            task_start_date = '171129';
        case '3358884-L'
            Threshold_Low = [-1.6,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6,-1.65,-1.6,-1.6,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [5,5,5,3,3,4,2,2,2,2,2,4];
            task_start_date = '171213';
        case '3358884-R'
            Threshold_Low = [-1.6,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6,-1.65,-1.6,-1.6,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [3,2,2,2,3,2,2,2,2,2,2,4];
            task_start_date = '171213';
        case '3358883-O'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.6,-1.65]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,6,3,2,3,2,4,2,3,2,2,4];
            task_start_date = '171214';
        case '3358883-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.6,-1.65]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [3,2,2,2,2,2,2,2,3,2,2,4];
            task_start_date = '171214';
        case '3373693-O'
            Threshold_Low = [-1.65,-1.6,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,2,2,2,2,2,3,2,2,2,2,4];
            task_start_date = '180111';
        case '3373693-L'
            Threshold_Low = [-1.65,-1.6,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.65,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [6,4,3,3,3,2,2,2,2,2,2,4];
            task_start_date = '180111';
        case '3373693-LR'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [4,4,2,2,2,2,2,2,2,1,2,7];
            task_start_date = '180117';
        case '3373693-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [3,3,2,3,2,4,2,2,2,2,2,4];
            task_start_date = '180117';   
        case '3420509-L'
            Threshold_Low = [-1.6,-1.65,-1.6,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [5,3,3,1,2,2,2,2,2,2,2];
            task_start_date = '180206';  
        case '3420509-R'
            Threshold_Low = [-1.6,-1.65,-1.6,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [5,5,3,3,3,2,2,2,2,2,2];
            task_start_date = '180206';   
        case '3358845-O'
            Threshold_Low = [-1.6,-1.65,-1.65,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [4,4,2,2,4,2,2,2,2,2,2];
            task_start_date = '180313';  
        case '3358845-LR'
            Threshold_Low = [-1.6,-1.65,-1.65,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [4,3,3,2,4,3,2,2,2,2,2];
            task_start_date = '180313';
        case '3438483-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,2,2,2,2,2,2,2,2,2];
            task_start_date = '180501'; 
        case '3491479-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [5,3,2,2,2,2,2,2,2,2,2];
            task_start_date = '180501'; 
        case '3491479-LR'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [5,3,3,2,2,2,2,2,2,2,3];
            task_start_date = '180503';
        case '3491479-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [5,4,2,2,2,2,2,2,2,2,2];
            task_start_date = '180508'; 
        case '3438500-O'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [4,5,2,2,2,2,2,2,2,2,2];
            task_start_date = '180508';
        case '3438500-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,4,2,2,2,2,2,2,2,2,2];
            task_start_date = '180508';
        case '3438521-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,3,2,3,2,2,2,2,2,2];
            task_start_date = '180530';
        case '3438521-LR'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [4,3,2,2,2,2,2,2,2,2,2];
            task_start_date = '180530';
        case '3438521-O'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,2,2,2,2,2,3,2,2,2];
            task_start_date = '180531';
        case '3438521-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [4,3,2,2,2,2,2,2,2,2,2];
            task_start_date = '180531';
        case '3453262-O'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,3,3,2,2,2,2,2,2,3,2];
            task_start_date = '180625';
        case '3453262-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,2,2,2,2,2,2,2,2,2];
            task_start_date = '180625';
        case '3438544-O'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,3,3,3,2,2,2,2,2,2,2];
            task_start_date = '180624';
        case '3438544-L'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,2,2,2,2,2,2,2,2,2];
            task_start_date = '180624';
        case '3438544-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [5,4,4,4,3,3,2,3,2,2,2];
            task_start_date = '180624';
        case '3495100-R'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [5,4,4,3,4,3,4,2,4,3,2];
            task_start_date = '180920';
        case '3547207-LR'
            Threshold_Low = [-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.65,-1.6,-1.6,-1.6,-1.6]; Threshold_High = Threshold_Low+0.2;
            ImBlock = [3,2,2,2,2,2,2,2,2,2,2];
            task_start_date = '181103';
        case '4383182-O'
            Threshold_Low = [-1.88,-1.88]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [3,2];
            task_start_date = '210904';
        case '4383182-L'
            Threshold_Low = [-1.88,-1.88]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [2,2];
            task_start_date = '210904';
        case '4383183-O'
            Threshold_Low = [-1.88,-1.88]; Threshold_High = Threshold_Low+0.25;
            ImBlock = [2,2];
            task_start_date = '210904';
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
    
    index = find(ismember(Date_list,task_start_date));
    Date_list = Date_list(index:end);
    
    % Load trial index
    load(['Z:\People\Chi\WFLP_IN\CheckNonconsTrials' filesep IN '_' Initial '_' Animal '_CheckNonconsTrials_Corrected.mat'], 'TrialInfo_xsg_all');
    TrialInfo_xsg_all = TrialInfo_xsg_all(index:end);
    %% Load Dispatcher and xsg files for current day
    disp([Initial '_' Animal ' ' num2str(Gap)]);

    % Initialization
    IndexInfo = cell(size(ImBlock))';
    IndexInfo_ITI = cell(size(ImBlock))';

    % Basically go throuh trial by trial
    for curr_Imaging_Day = 1:length(ImBlock)
            
        curr_day = Imaging_Day(curr_Imaging_Day);
        ImBlockNum = ImBlock(curr_Imaging_Day);
        TrialInfo_xsg = TrialInfo_xsg_all{curr_day};
        IndexInfo{curr_Imaging_Day,1} = cell(ImBlockNum,1);
        IndexInfo_ITI{curr_Imaging_Day,1} = cell(ImBlockNum,1);
        
        Date = Date_list{curr_day}; disp([num2str(curr_day) ': ' Date]);
        
        %% Dispatcher, should only have one filename
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
        TrialNum_Total = length(TrialInfo_Dispatcher); 
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
            continue
        end
        Xsg_filename_list = sort(Xsg_filename_list); % from 1st to end


        %% Load xsg file sequentially

        for curr_ImBlock = 1:ImBlockNum
            if curr_Imaging_Day <= 11
                curr_xsg = 2*curr_ImBlock-1;
            else
                curr_xsg = curr_ImBlock;
            end
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
            
            consec_trials = diff([0; TrialInfo_xsg{curr_xsg}(:,2)]);
            if sum(consec_trials(2:end)~=1) > 1
                disp('Non-consecutive trials');
            end
            
            % Ignore if no trials exist
            if isempty(TrialInfo_xsg{curr_xsg})
                disp(['No Bitecode in ' Xsg_filename_list{curr_xsg} '.'])
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
            [Lever_active, Lever_force_resample Lever_force_smooth ...
                Lever_velocity_envelope_smooth] ...
                = CR_parseLeverMovement_Updated(Lever_trace_xsg, Gap, Leeway); % 1K Hz
            Lever_active_switch = diff([0;Lever_active]); % 1K Hz
            
            % Get movment start and stop from lever trace
            [Licking_active, Licking_bout, Licking_bout_switch] = CR_GetLickingState(Licking_xsg,3); % 1K Hz
         
            for curr_trial = 1:size(TrialInfo_xsg{curr_xsg},1)
                % Look at all trials
                curr_trial_index = TrialInfo_xsg{curr_xsg}(curr_trial,2);
                if curr_trial_index > TrialNum_Total
                    continue
                end
                Bitcode_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.bitcode(1);
                % Cue time
                Cue_Start_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.cue(1);
                Cue_End_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.cue(2);
                Cue_Start_xsg = Cue_Start_dispatcher - Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                Cue_End_xsg = Cue_End_dispatcher - Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                % Reward time
                if ismember(curr_trial_index, Rewarded_trial)
                    Reward_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.reward(1);
                    Reward_xsg = Reward_dispatcher - Bitcode_dispatcher + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                else
                    Reward_xsg = nan;
                end
                % ITI
                ITI_start_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.iti(1);
                ITI_end_dispatcher = TrialInfo_Dispatcher{curr_trial_index,1}.states.iti(2);
                ITI_start_xsg = ITI_start_dispatcher - Bitcode_dispatcher ...
                    + TrialInfo_xsg{curr_xsg}(curr_trial,1);
                ITI_end_xsg = ITI_end_dispatcher - Bitcode_dispatcher ...
                    + TrialInfo_xsg{curr_xsg}(curr_trial,1);                     
                % End Control
                End_xsg = max(Cue_End_xsg, Reward_xsg);
                % Refine, within 300 sec
                if Cue_Start_xsg>300 || End_xsg>300
                    continue
                end
                
                Temp_FrameIndex_Info = nan(1,13);
                temp_movonset = [];
                
                Temp_FrameIndex_Info(1,1) = curr_trial_index; % Column 1: Trial index;
                Temp_FrameIndex_Info(1,2) = Cue_Start_xsg; % Column 2: Cue Onset;
                Temp_FrameIndex_Info(1,13) = CatchInfo_Dispatcher(curr_trial_index);
                if sum(Lever_active((round(Cue_Start_xsg*1000)-100):round(Cue_Start_xsg*1000)))==0
                    Temp_FrameIndex_Info(1,3) = 1; % Column 3: Cued or not;
                    temp_movonset = find(Lever_active_switch(round(Cue_Start_xsg*1000):round(End_xsg*1000))==1,1);
                    if ~isempty(temp_movonset)
                        Temp_FrameIndex_Info(1,4) = 1; % Column 4: Responded or not;
                        Temp_FrameIndex_Info(1,5) = temp_movonset./1000; % Column 5: Reaction Time;
                    end
                end
                Temp_FrameIndex_Info(1,6) = Reward_xsg; % Column 6: Reward Time;
                if ~isnan(Reward_xsg)
                    curr_sample_range = [];
                    temp_lever_trace = [];
                    temp_velocity = [];
                    temp_crossing_higher = [];
                    temp_crossing_lower = [];
                    if isempty(find(Lever_active_switch(1:round(Reward_xsg*1000))==1,1,'last'))
                        Temp_FrameIndex_Info(1,6) = nan;
                        continue
                    elseif size(IndexInfo{curr_Imaging_Day,1}{curr_ImBlock,1},1)>0 && find(Lever_active_switch(1:round(Reward_xsg*1000))==1,1,'last') == (IndexInfo{curr_Imaging_Day,1}{curr_ImBlock,1}(end,7)*1000)
                        disp([Animal ' ' num2str(curr_Imaging_Day) ' ' num2str(curr_ImBlock)])
                        Temp_FrameIndex_Info(1,6) = nan;
                        continue
                    end
                    Temp_FrameIndex_Info(1,7) = find(Lever_active_switch(1:round(Reward_xsg*1000))==1,1,'last')/1000; % Column 7: Reward MovOnst Time;
                    % Column 8: Movement Duration
                    temp_mov_duration = find(Lever_active_switch(round(Temp_FrameIndex_Info(1,7)*1000):end)==-1,1,'first')/1000;
                    if isempty(temp_mov_duration);
                        temp_mov_duration = 300 - Temp_FrameIndex_Info(1,7);
                    end
                    Temp_FrameIndex_Info(1,8) = temp_mov_duration;
                    % Crossing Low and High
                    curr_sample_range = round(Temp_FrameIndex_Info(1,7)*1000)+[0:round(Temp_FrameIndex_Info(1,8)*1000)];
                    temp_lever_trace = Lever_force_resample(curr_sample_range);
                    temp_velocity = diff([0;temp_lever_trace]);                   
                    temp_crossing_higher = find(temp_lever_trace <= -Threshold_High(curr_Imaging_Day)-2.5132);
                    if ~isempty(temp_crossing_higher)
                        for curr_crossing_higher = 1:length(temp_crossing_higher)
                            if temp_velocity(temp_crossing_higher(curr_crossing_higher))<=0
                                temp_crossing_lower = find(temp_lever_trace(1:temp_crossing_higher(curr_crossing_higher)) >= -Threshold_Low(curr_Imaging_Day)-2.5132,1,'last');
                                if ~isempty(temp_crossing_lower)
                                    % Column 9: Time crossing lower threshold
                                    Temp_FrameIndex_Info(1,9) = (temp_crossing_lower-1)./1000+Temp_FrameIndex_Info(1,7);
                                    temp_crossing_higher = temp_crossing_higher(curr_crossing_higher);
                                    % Time crossing higher threshold
                                    Temp_FrameIndex_Info(1,10) = (temp_crossing_higher-1)./1000+Temp_FrameIndex_Info(1,7);
                                    break
                                end
                            end
                        end
                    end
                    % Licking information
                    temp_lick_onset = [];
                    temp_lick_duration = [];
                    End_xsg = min(ITI_end_xsg, 300);
                    if ~isempty(Licking_bout_switch)
                        temp_lick_onset = find(Licking_bout_switch(round(Temp_FrameIndex_Info(1,7)*1000):round(End_xsg*1000))==1,1,'first');
                        if ~isempty(temp_lick_onset)
                            Temp_FrameIndex_Info(1,11) = (temp_lick_onset-1)/1000;
                            temp_lick_duration = find(Licking_bout_switch(round(Temp_FrameIndex_Info(1,7)*1000)+temp_lick_onset-1:round(End_xsg*1000))==-1,1,'first')-1;
                            if ~isempty(temp_lick_duration)
                                Temp_FrameIndex_Info(1,12) = temp_lick_duration/1000;
                            else
                                Temp_FrameIndex_Info(1,12) = inf(1);
                            end
                        elseif isempty(temp_lick_onset) && Licking_bout(round(Temp_FrameIndex_Info(1,7)*1000)) == 1
                            % lick starts within 500ms before moveonset
                            temp_lick_onset = find(Licking_bout_switch(round(Temp_FrameIndex_Info(1,7)*1000-500):round(End_xsg*1000))==1,1,'first');
                            if ~isempty(temp_lick_onset)
                                Temp_FrameIndex_Info(1,11) = (temp_lick_onset-1)/1000-0.5;
                                temp_lick_duration = find(Licking_bout_switch(round(Temp_FrameIndex_Info(1,7)*1000)+temp_lick_onset-1:round(End_xsg*1000))==-1,1,'first')-1;
                                if ~isempty(temp_lick_duration)
                                    Temp_FrameIndex_Info(1,12) = temp_lick_duration/1000;
                                else
                                    Temp_FrameIndex_Info(1,12) = inf(1);
                                end
                            else
                                Temp_FrameIndex_Info(1,11) = inf(1);
                                Temp_FrameIndex_Info(1,12) = inf(1);
                            end
                        end
                        else
                        Temp_FrameIndex_Info(1,11) = inf(1);
                        Temp_FrameIndex_Info(1,12) = inf(1);
                    end
                end
                
                IndexInfo{curr_Imaging_Day,1}{curr_ImBlock,1} = vertcat(IndexInfo{curr_Imaging_Day,1}{curr_ImBlock,1},Temp_FrameIndex_Info);
                
                % ITI movement
                temp_movonset = [];
                if ITI_start_xsg >= 300
                    continue
                end
                End_xsg = min(ITI_end_xsg, 300);
                temp_movonset = find(Lever_active_switch(round(ITI_start_xsg*1000):round(End_xsg*1000))==1);
                if ~isempty(temp_movonset)
                    for curr_mov = 1:length(temp_movonset)                                                
                        curr_movonset = temp_movonset(curr_mov);
                        temp_mov_duration = find(Lever_active_switch(round(ITI_start_xsg*1000)+curr_movonset-1:round(End_xsg*1000))==-1,1,'first');
                        if isempty(temp_mov_duration)
                            continue
                        end
                        curr_sample_range = [];
                        temp_lever_trace = [];
                        temp_velocity = [];
                        temp_crossing_higher = [];
                        temp_crossing_lower = [];
                        curr_sample_range = round(ITI_start_xsg*1000)+curr_movonset-1+[0:temp_mov_duration];
                        temp_lever_trace = Lever_force_resample(curr_sample_range);
                        temp_velocity = diff([0;temp_lever_trace]);                   
                        temp_crossing_higher = find(temp_lever_trace <= -Threshold_High(curr_Imaging_Day)-2.5132);
                        if ~isempty(temp_crossing_higher)
                            for curr_crossing_higher = 1:length(temp_crossing_higher)
                                if temp_velocity(temp_crossing_higher(curr_crossing_higher))<=0
                                    temp_crossing_lower = find(temp_lever_trace(1:temp_crossing_higher(curr_crossing_higher)) >= -Threshold_Low(curr_Imaging_Day)-2.5132,1,'last');
                                    if ~isempty(temp_crossing_lower)
                                        % Column 9: Time crossing lower threshold
                                        Temp_FrameIndex_Info_ITI = nan(1,7);
                                        Temp_FrameIndex_Info_ITI(1,1) = curr_trial_index;
                                        Temp_FrameIndex_Info_ITI(1,2) = ITI_start_xsg+(curr_movonset-1)/1000;
                                        Temp_FrameIndex_Info_ITI(1,3) = temp_mov_duration./1000;
                                        Temp_FrameIndex_Info_ITI(1,4) = (temp_crossing_lower-1)/1000+Temp_FrameIndex_Info_ITI(1,2);
                                        Temp_FrameIndex_Info_ITI(1,5) = (temp_crossing_higher(curr_crossing_higher)-1)/1000+Temp_FrameIndex_Info_ITI(1,2);
                                        % Licking
                                        temp_lick_onset = [];
                                        temp_lick_duration = [];
                                        if curr_mov < length(temp_movonset)
                                           End_control = (temp_movonset(curr_mov+1)-1)/1000+ITI_start_xsg;
                                        else
                                           End_control = End_xsg;
                                        end
                                        if ~isempty(Licking_bout_switch)
                                            temp_lick_onset = find(Licking_bout_switch(round(Temp_FrameIndex_Info_ITI(1,2)*1000):round(End_control*1000))==1,1,'first');
                                            if ~isempty(temp_lick_onset)
                                                Temp_FrameIndex_Info_ITI(1,6) = temp_lick_onset/1000;
                                                temp_lick_duration = find(Licking_bout_switch(round(Temp_FrameIndex_Info_ITI(1,2)*1000)+temp_lick_onset-1:round(End_xsg*1000))==-1,1,'first');
                                                if ~isempty(temp_lick_duration)
                                                    Temp_FrameIndex_Info_ITI(1,7) = temp_lick_duration/1000;
                                                end
                                            elseif isempty(temp_lick_onset) && Licking_bout(round(Temp_FrameIndex_Info_ITI(1,2)*1000)) == 1
                                                Temp_FrameIndex_Info_ITI(1,6) = inf(1);
                                                Temp_FrameIndex_Info_ITI(1,7) = inf(1);
                                            end
                                        else
                                            Temp_FrameIndex_Info_ITI(1,6) = inf(1);
                                            Temp_FrameIndex_Info_ITI(1,7) = inf(1);
                                        end
                                        IndexInfo_ITI{curr_Imaging_Day,1}{curr_ImBlock,1} = vertcat(IndexInfo_ITI{curr_Imaging_Day,1}{curr_ImBlock,1},Temp_FrameIndex_Info_ITI);
                                        break
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            disp(['Finish ' Xsg_filename_list{curr_xsg}]);

        end
        %
        disp('Finish all xsg files');
        disp(['Finish ' Date]);
    end

    %% Save

    % User interaction
    disp('Saving...');

    % Save
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'MovAnalysis']);

    curr_filename_lever = [Initial '_' Animal '_FrameIndex_Info'];
    save(curr_filename_lever, 'IndexInfo','IndexInfo_ITI','Imaging_Day','-v7.3');

    % User interaction
    disp('Saving done');

end

disp('Finish All animals');
    

