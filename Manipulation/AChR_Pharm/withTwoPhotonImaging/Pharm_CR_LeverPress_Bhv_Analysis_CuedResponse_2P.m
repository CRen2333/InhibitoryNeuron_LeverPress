%% ***** Formal Coding V2.0 *****

% Frame 13 is movement onset
clear;
close all;
clc;

IN = 'VIP';
% Animals = {'CR_3702608-O','CR_3702608-LR'};
% Animals = {'CR_3658844-L','CR_3658844-R'};
% Animals = {'CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O'};
% Animals = {'CR_3887044-O','CR_3887044-L'};
% Animals = {'CR_3786160-LR'};
% Animals = {'CR_3887043-O','CR_3887043-L','CR_3887041-O'};
% SessionType = 'Nai_Ant';
% Animals = {'CR_3672031-O'};
% SessionType = 'Exp_Ago_Ant';
% Animals = {'CR_3702608-LR'};
% SessionType = 'Exp_Ant';
% SessionType = 'Exp_Ago';
% Animals = {'WL_3526642-R'};
% SessionType = 'Nai_nAnt';
% Animals = {'CR_4259757-R','CR_4302983-LL','CR_4302983-R'};
% Animals = {'CR_4302984-O','CR_4302984-R','CR_4302984-LR'};
% Animals = {'CR_4302984-R'};
SessionType = 'Nai_mAnt';
% Animals = {'CR_4303048-L','CR_4303048-R','CR_4303048-LR'};
Animals = {'CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'};
% % IN = 'SOM';
% Animals = {'CR_3672035-L','CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O'};
% Animals = {'CR_3886982-O','CR_3886982-L'}; 
% Animals = {'CR_3886982-R','CR_3886982-LR'}; 
% Animals = {'CR_3887040-O','CR_3887040-L'};
% SessionType = 'Nai_Ant';
% Animals = {'WL_3547273-LR','WL_3547272-O','WL_3547272-L','WL_3526578-O','WL_3526578-R'}; 
% Animals = {'WL_3547272-O'};
% SessionType = 'Exp_Ago';
% Animals = {'CR_3672035-L','CR_3672035-R'};
% Animals = {'CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O'};
% SessionType = 'Exp_Ago_Ant';

Gap = round(3*33.33)*5;
MovDur = 200;
Leeway = 150;

for curr_animal = 1:length(Animals)

    %% Clean
    clearvars -except Initial Animals IN curr_animal Gap MovDur Leeway SessionType;
    clc;
    Animal = Animals{curr_animal};
    disp(Animal);
    
    % Set sample range
    Sample_range_xsg = [-5000:40000]; % 1 ms = 10 frames in xsg;
    Sample_range_resample = [-500:4000];
    Sample_range_corr = [50:50+MovDur]; % 10 ms step
    Lever_downsamp = 10;

    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep SessionType];
    cd(General_path);
    mkdir('MovAnalysis');

    %% Load files    
    % Get Dates
    Xsg_path = [General_path filesep 'LeverTrace'];
    Dates = dir(Xsg_path);
    Dates = {Dates.name};
    Dates = sort(Dates);
    Dates = Dates(3:end);
    % For Bpod
    Bpod_path = [General_path filesep 'Bpod'];
    Bpod_list = dir(Bpod_path);
    Bpod_list = {Bpod_list.name};
    Bpod_list = Bpod_list(3:end);

    % Initialization
    CuedMov_SingleAnimal = cell(size(Dates));
    CuedRewardedFrameIndex = cell(size(Dates));
    MovActiveFrame = cell(size(Dates));
    MovActiveFrame_PreExtended = cell(size(Dates));

    % Basically go throuh trial by trial
    for curr_day = 1:length(Dates)
        % User interaction
        Date = Dates{curr_day};
        disp(['Running ' Animal ' on ' Date]);
        clear Frame_times Lever_active Lever_active_switch Lever_force_resample TrialInfo_xsg
        clear TrialInfo_Bpod

        Curr_Bpod_filname = Bpod_list(cellfun(@(x) ~isempty(strfind(x, Date(3:end))), Bpod_list));
        if isempty(Curr_Bpod_filname)
            warning(['No Dispatcher file saved for ' Animal ' on ' Date '.']);
            continue
        end
        % Load extracted xsg channels
        load([Xsg_path filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times','Lever_active','Lever_force_resample','TrialInfo_xsg','-mat');
        for curr_session = 1:length(Curr_Bpod_filname)
            warning off
            load([Bpod_path filesep Curr_Bpod_filname{curr_session}]);
            warning on
        
            % Get general trial information that will be used for all xsg files
            % Extract trial information from Bpod
            TrialInfo_Bpod = SessionData.RawEvents.Trial';
            TrialStart_Bpod = SessionData.TrialStartTimestamp';
            clear SessionData;
            % Total trial num
            TrialNum_Total = min(length(TrialInfo_Bpod),size(TrialInfo_xsg{curr_session},1));
            temp_Bpod = TrialStart_Bpod; % column 1: bitcode time
            for curr_trial = 1:TrialNum_Total
                % column 2: ITI start
                temp_Bpod(curr_trial,2) = temp_Bpod(curr_trial,1) + TrialInfo_Bpod{curr_trial,1}.States.ITI(1);
                % column 3: cue start
                temp_Bpod(curr_trial,3) = temp_Bpod(curr_trial,1) + TrialInfo_Bpod{curr_trial,1}.States.Cue(1);
                % column 4: cue end
                temp_Bpod(curr_trial,4) = temp_Bpod(curr_trial,1) + TrialInfo_Bpod{curr_trial,1}.States.CueTup(1);
                % column 5: reward
                temp_Bpod(curr_trial,5) = temp_Bpod(curr_trial,1) + TrialInfo_Bpod{curr_trial,1}.States.Reward(1);
            end
            Time_offset = temp_Bpod(1:TrialNum_Total,1) - TrialInfo_xsg{curr_session}(1:TrialNum_Total,1);            
            temp_Bpod = temp_Bpod(1:TrialNum_Total,:) - repmat(Time_offset,1,5);
            clear TrialInfo_Bpod
            TrialInfo_Bpod{curr_session} = temp_Bpod;
            clear Time_offset temp_Bpod TrialStart_Bpod
        
            % Get rewarded trials
            Rewarded_trial = find(~isnan(TrialInfo_Bpod{curr_session}(:,5)));

            % Initialization for each day
            CuedMov.Cued_MovOnset_Info{curr_session} = [];
            CuedMov.LeverTrace{curr_session} = [];
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
            CuedMov.CR = nan;

            % Get movment start and stop from lever trace
            Lever_active_switch{curr_session} = diff([0;Lever_active{curr_session}]); % 1K Hz

            MovOnset_curr_session = [];
            LeverTrace_curr_session = [];

            % Get frame onset information
            if ~isempty(Frame_times{curr_session})
                TimeUp = Frame_times{curr_session}(end);
                if TimeUp > length(Lever_active{curr_session})
                    disp('Imaging exceeds behavior!')
                    TimeUp = length(Lever_active{curr_session});
                    Frame_times{curr_session}(Frame_times{curr_session}>TimeUp) = [];
                end
            else
                TimeUp = length(Lever_active{curr_session});
            end

            for curr_trial = 1:TrialNum_Total

                % Cue time index in lever trace
                Cue_Start_xsg = round(TrialInfo_Bpod{curr_session}(curr_trial,3)*1000);
                Cue_End_xsg = round(TrialInfo_Bpod{curr_session}(curr_trial,4)*1000);
                Reward_xsg = round(TrialInfo_Bpod{curr_session}(curr_trial,5)*1000);
                if isnan(Cue_End_xsg)
                    End_xsg = Reward_xsg;
                else
                    End_xsg = Cue_End_xsg;
                end
                % Refine
                if End_xsg > TimeUp*1000
                    continue
                end

                Temp_MovOnset_info = [];
                Temp_MovTrace_curr_trial = [];
                MovOnset_curr_trial = [];
                MovOnset_curr_trial = find(Lever_active_switch{curr_session}(round(Cue_Start_xsg):round(End_xsg))==1,1);
                % Ignore if no movement during the cue
                if ~isempty(MovOnset_curr_trial)
                    disp(num2str(curr_trial));
                    MovOnset_curr_trial = find(Lever_active_switch{curr_session}(round(Cue_Start_xsg):round(End_xsg))==1) + round(Cue_Start_xsg) - 1;
                    Rewarded_MovOnset = nan;
                    % Rewarded or not
                    if ~isnan(Reward_xsg)
                        Temp_Rewarded_MovOnset = find(Lever_active_switch{curr_session}(1:round(Reward_xsg))==1,1,'last'); % index
                        % Cued or not
                        if ~isempty(Temp_Rewarded_MovOnset)
                            if Temp_Rewarded_MovOnset >= Cue_Start_xsg && sum(Lever_active{curr_session}((round(Cue_Start_xsg)-100):round(Cue_Start_xsg)))==0 && (Temp_Rewarded_MovOnset+Sample_range_resample(end))/1000 <= TimeUp % Collect trace
                                Rewarded_MovOnset = Temp_Rewarded_MovOnset;
                            end
                        end
                    end

                    for curr_movonset = 1:length(MovOnset_curr_trial)
                        % ignore if exceeds
                        if (MovOnset_curr_trial(curr_movonset)+Sample_range_resample(1))/1000 <0 || (MovOnset_curr_trial(curr_movonset)+Sample_range_resample(end))/1000 >TimeUp
                            continue
                        end
                        Temp_MovOnset_info(curr_movonset,1) = curr_trial;
                        Temp_MovOnset_info(curr_movonset,2) = MovOnset_curr_trial(curr_movonset)/1000; % sec
                        % Get movement duration
                        Temp_MovOnset_info(curr_movonset,3) = find(Lever_active_switch{curr_session}(MovOnset_curr_trial(curr_movonset):end)==-1,1,'first')/1000;
                        if isnan(Temp_MovOnset_info(curr_movonset,3))
                            Temp_MovOnset_info(curr_movonset,3) = TimeUp - Temp_MovOnset_info(curr_movonset,2);
                        end
                        Temp_MovOnset_info(curr_movonset,4) = curr_movonset;
                        % Cued_rewarded movement or not
                        if ~isnan(Rewarded_MovOnset)
                            if MovOnset_curr_trial(curr_movonset) == Rewarded_MovOnset
                                Temp_MovOnset_info(curr_movonset,5) = true;
                                % 500ms before cued rewarded movementonset is quiet or not
                                if sum(Lever_active{curr_session}((round(Rewarded_MovOnset)-500):round(Rewarded_MovOnset-1)))==0
                                    Temp_MovOnset_info(curr_movonset,8) = true;
                                else
                                    Temp_MovOnset_info(curr_movonset,8) = false;
                                end
                                % 1s before cued rewarded movementonset is quiet or not
                                if round(Rewarded_MovOnset)-1000 >= 0
                                    if sum(Lever_active{curr_session}((round(Rewarded_MovOnset)-1000):round(Rewarded_MovOnset-1)))==0
                                        Temp_MovOnset_info(curr_movonset,9) = true;
                                    else
                                        Temp_MovOnset_info(curr_movonset,9) = false;
                                    end
                                else
                                    Temp_MovOnset_info(curr_movonset,9) = false;
                                end
                            else
                                Temp_MovOnset_info(curr_movonset,5) = false;
                                Temp_MovOnset_info(curr_movonset,8) = false;
                                Temp_MovOnset_info(curr_movonset,9) = false;
                            end
                        else
                            Temp_MovOnset_info(curr_movonset,8) = false;
                            Temp_MovOnset_info(curr_movonset,9) = false;
                        end
                        Temp_MovOnset_info(curr_movonset,6) = Cue_Start_xsg/1000;
                        if isnan(Reward_xsg)
                            Temp_MovOnset_info(curr_movonset,7) = nan;
                        else
                            Temp_MovOnset_info(curr_movonset,7) = Reward_xsg/1000;
                        end
                        Temp_MovTrace_curr_trial = horzcat(Temp_MovTrace_curr_trial,Lever_force_resample{curr_session}(round(MovOnset_curr_trial(curr_movonset)+Sample_range_resample)));
                    end
                    MovOnset_curr_session = vertcat(MovOnset_curr_session,Temp_MovOnset_info);
                    LeverTrace_curr_session = horzcat(LeverTrace_curr_session,Temp_MovTrace_curr_trial);
                end
            end
            CuedMov.Cued_MovOnset_Info{curr_session} = MovOnset_curr_session;
            CuedMov.LeverTrace{curr_session} = LeverTrace_curr_session;
            CuedMov.Bpod{curr_session} = TrialInfo_Bpod{curr_session};
            CuedMov.CR = length(Rewarded_trial)/TrialNum_Total;
        end
        % All
        CuedMov.Cued_MovOnset_Info_All = cell2mat(CuedMov.Cued_MovOnset_Info');
        CuedMov.LeverTrace_All = cell2mat(CuedMov.LeverTrace);
        if isempty(CuedMov.LeverTrace_All)
            CuedMov.LeverTrace_All_downsample = [];
        else
            Temp_trace = resample(CuedMov.LeverTrace_All,1,Lever_downsamp);
            CuedMov.LeverTrace_All_downsample = Temp_trace(Sample_range_corr,:);
        end
        % First movement
        if isempty(CuedMov.LeverTrace_All)
            CuedMov.LeverTrace_First = [];        
            CuedMov.LeverTrace_First_downsample = [];
        else
            First_index = find(CuedMov.Cued_MovOnset_Info_All(:,4)==1);
            CuedMov.LeverTrace_First = CuedMov.LeverTrace_All(:,First_index);        
            Temp_trace = resample(CuedMov.LeverTrace_First,1,Lever_downsamp);
            CuedMov.LeverTrace_First_downsample = Temp_trace(Sample_range_corr,:);
        end
        % Reward movement
        
        if isempty(CuedMov.LeverTrace_All)
            CuedMov.LeverTrace_First_downsample = [];
             CuedMov.LeverTrace_Reward = [];
        else
            Reward_index = logical(CuedMov.Cued_MovOnset_Info_All(:,5));
            CuedMov.LeverTrace_Reward = CuedMov.LeverTrace_All(:,Reward_index);
            Temp_trace = resample(CuedMov.LeverTrace_Reward,1,Lever_downsamp);
            % CuedMov.LeverTrace_Reward_downsample = Temp_trace(Sample_range_corr_R,:);
            CuedMov.LeverTrace_Reward_downsample = Temp_trace(Sample_range_corr,:);
        end
        
        if ~isempty(CuedMov.Cued_MovOnset_Info_All)
            CuedMov.Median_Movduration = median(CuedMov.Cued_MovOnset_Info_All(:,3));
            CuedMov.Median_Movduration_Reward = median(CuedMov.Cued_MovOnset_Info_All(Reward_index,3));       
            CuedMov.Median_Movduration_First = median(CuedMov.Cued_MovOnset_Info_All(First_index,3));
        end
        CuedMov_SingleAnimal{curr_day} = CuedMov;
        
        for curr_session = 1:length(Frame_times)
            if isempty(CuedMov.Cued_MovOnset_Info{curr_session})
                CuedRewardedFrameIndex{curr_day}{curr_session} = [];
                CuedRewardedFrameIndex_z1{curr_day}{curr_session} = [];
                CuedRewardedFrameIndex_z2{curr_day}{curr_session} = [];
                continue
            end
            Index = CuedMov.Cued_MovOnset_Info{curr_session}(:,5);
            CuedRewardedMovOnset = CuedMov.Cued_MovOnset_Info{curr_session}(logical(Index),2);
            % Super time consuming! Try this!
            MovActiveFrame{curr_day}{curr_session} = Lever_active{curr_session}(round(Frame_times{curr_session}*1000));
            % Get Frame index
            if ~isempty(Frame_times{curr_session})
                for curr_trial = 1:length(CuedRewardedMovOnset)
                    curr_onset = CuedRewardedMovOnset(curr_trial);
                    [C,I_onset] = min(abs(Frame_times{curr_session}-curr_onset));
                    CuedRewardedFrameIndex{curr_day}{curr_session}(:,curr_trial) = [I_onset-14:I_onset+57];
                    Frame_times_z1{curr_session} = Frame_times{curr_session}(1:2:length(Frame_times{curr_session}));
                    [C,I_onset_z1] = min(abs(Frame_times_z1{curr_session}-curr_onset));
                    CuedRewardedFrameIndex_z1{curr_day}{curr_session}(:,curr_trial) = [I_onset_z1-7:I_onset_z1+28];
                    Frame_times_z2{curr_session} = Frame_times{curr_session}(2:2:length(Frame_times{curr_session}));
                    [C,I_onset_z2] = min(abs(Frame_times_z2{curr_session}-curr_onset));
                    CuedRewardedFrameIndex_z2{curr_day}{curr_session}(:,curr_trial) = [I_onset_z2-7:I_onset_z2+28];
                end
                clear Frame_times_z1 Frame_times_z2
            else
                CuedRewardedFrameIndex{curr_day}{curr_session} = [];
                CuedRewardedFrameIndex_z1{curr_day}{curr_session} = [];
                CuedRewardedFrameIndex_z2{curr_day}{curr_session} = [];
            end

            % Get index for 500 and 1000ms quiet
            if isempty(CuedMov.Cued_MovOnset_Info_All)
                Index_1000{curr_day}{curr_session} = [];
                Index_500{curr_day}{curr_session} = [];
            else
                Index = logical(CuedMov.Cued_MovOnset_Info_All(:,5));
                Index_1000{curr_day}{curr_session} = logical(CuedMov.Cued_MovOnset_Info_All(Index,9));
                Index_500{curr_day}{curr_session} = logical(CuedMov.Cued_MovOnset_Info_All(Index,8));
            end
        end

        disp('Finish all sessions');
        disp(['Finish ' Date]);
    end

    %% Calculate Correlation Matrix
    disp('Calculating correlation...')
    % Trail_Trial_Correlation Matrix
    % All traces
    Trial_Trial_Corr_All = nan(length(Dates));
    for curr_day_1 = 1:length(Dates)
        if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample)
            for curr_day_2  = 1:length(Dates)
                if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_All_downsample) 
                    if curr_day_1 == curr_day_2
                        Trial_Trial_Corr_All(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample);
                    else
                        temp_matrix = [];
                        temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_All_downsample, ...
                            CuedMov_SingleAnimal{curr_day_2}.LeverTrace_All_downsample);
                        Trial_Trial_Corr_All(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(temp_matrix);
                    end
                end
            end
        end
    end

    % First traces
    Trial_Trial_Corr_First = nan(length(Dates));
    for curr_day_1 = 1:length(Dates)
        if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample)
            for curr_day_2  = 1:length(Dates)
                if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_First_downsample) 
                    if curr_day_1 == curr_day_2
                        Trial_Trial_Corr_First(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample);
                    else
                        temp_matrix = [];
                        temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_First_downsample, ...
                            CuedMov_SingleAnimal{curr_day_2}.LeverTrace_First_downsample);
                        Trial_Trial_Corr_First(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(temp_matrix);
                    end
                end
            end
        end
    end

    % Reward traces
    Trial_Trial_Corr_Reward = nan(length(Dates));
    for curr_day_1 = 1:length(Dates)
        if ~isempty(CuedMov_SingleAnimal{curr_day_1}) && ~isempty(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample)
            for curr_day_2  = 1:length(Dates)
                if ~isempty(CuedMov_SingleAnimal{curr_day_2}) && ~isempty(CuedMov_SingleAnimal{curr_day_2}.LeverTrace_Reward_downsample) 
                    if curr_day_1 == curr_day_2
                        Trial_Trial_Corr_Reward(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample);
                    else
                        temp_matrix = [];
                        temp_matrix = horzcat(CuedMov_SingleAnimal{curr_day_1}.LeverTrace_Reward_downsample, ...
                            CuedMov_SingleAnimal{curr_day_2}.LeverTrace_Reward_downsample);
                        Trial_Trial_Corr_Reward(curr_day_1,curr_day_2) = CR_Get_Median_Corrcoef(temp_matrix);
                    end
                end
            end
        end
    end

    %% Save

    % User interaction
    disp('Saving...');

    % Save
    cd(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep SessionType filesep 'MovAnalysis']);
    curr_filename = [Animal '_CuedMov_All'];
    save(curr_filename, 'CuedMov_SingleAnimal','Trial_Trial_Corr_All','Trial_Trial_Corr_First','Trial_Trial_Corr_Reward','-v7.3');
    curr_filename = [Animal '_CuedRewardedMov_FrameIndex'];
    save(curr_filename, 'CuedRewardedFrameIndex','CuedRewardedFrameIndex_z1','CuedRewardedFrameIndex_z2','MovActiveFrame','-v7.3');
    curr_filename = [Animal '_goodTrialIndex'];
    save(curr_filename, 'Index_1000','Index_500','-v7.3');
    
    % User interaction
    disp('Saving done');
    disp(['Finish Gap ' num2str(Gap) 'ms.']);
end
disp('Finish All animals');




