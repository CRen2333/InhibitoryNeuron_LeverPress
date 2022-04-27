%% Collect_behavior_2P
clearvars -except VIPSOM_2P;
close all;
clc;

IN = 'VIP_SOM';
Animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R',...
    'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
Stage = 'Nai';

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage '\MovAnalysis\' Animal '_CuedMov_All.mat']);
    temp_matrix = nan(2,2);
    temp_matrix(1:size(Trial_Trial_Corr_Reward,1),1:size(Trial_Trial_Corr_Reward,1)) = Trial_Trial_Corr_Reward;
    LeverCorr_Reward_Matrix(:,:,curr_animal) = temp_matrix;
    RwdMov_Duration(curr_animal,[1:2]) = nan;
    CR(curr_animal,[1:2]) = nan;
    RT(curr_animal,[1:2]) = nan;
    RwdMVM_Rwd(curr_animal,[1:2]) = nan;
    C_CRM(curr_animal,[1:2]) = nan;
    TrialNum(curr_animal,[1:2]) = nan;
    
    for curr_day = 1:length(CuedMov_SingleAnimal)
        if isempty(CuedMov_SingleAnimal{curr_day})
            continue
        end
        RwdMov_Duration(curr_animal,curr_day) = CuedMov_SingleAnimal{curr_day}.Median_Movduration_Reward;
        temp_bpod = cell2mat(CuedMov_SingleAnimal{curr_day}.Bpod');
        CR(curr_animal,curr_day) = sum(~isnan(temp_bpod(:,5)))./size(temp_bpod,1);            
        temp_info = CuedMov_SingleAnimal{curr_day}.Cued_MovOnset_Info_All;
        first_index = temp_info(:,4)==1;
        reward_index = ~isnan(temp_info(:,7));
        rewared_mvm_index = temp_info(:,5)==1;
        RT(curr_animal,curr_day) = nanmedian(temp_info(first_index,2)-temp_info(first_index,6));
        RwdMVM_Rwd(curr_animal,curr_day) = nanmedian(temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2));
        C_CRM(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,2)-temp_info(rewared_mvm_index,6));
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
        temp_RwdMVM_Rwd = temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2);
        TrialNum(curr_animal,curr_day) = size(temp_trace,2);
        temp_RwdMVM_Rwd = round(temp_RwdMVM_Rwd*100);
        clear temp_RwdMVM_Rwd

    end
    clear clear temp_trace Trial_Trial_Corr_Reward CuedMov_SingleAnimal
end

for ii = 1:length(Animals)
    VIPSOM_2P.LeverCorr_Reward_within(ii,:) = diag(LeverCorr_Reward_Matrix(:,:,ii));
end

VIPSOM_2P.LeverCorr_Reward_Matrix = LeverCorr_Reward_Matrix;
VIPSOM_2P.RwdMov_Duration = RwdMov_Duration;
VIPSOM_2P.CR = CR;
VIPSOM_2P.RT = RT;
VIPSOM_2P.RwdMVM_Rwd = RwdMVM_Rwd;
VIPSOM_2P.C_CRM = C_CRM;
VIPSOM_2P.TrialNum = TrialNum;

% Check basic movement measurements
hM4Di_animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R'};
mCherry_animals = {'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
hM4Di_index = ismember(Animals,hM4Di_animals);
mCherry_index = ismember(Animals,mCherry_animals);

for curr_animal = 1:length(hM4Di_animals)
    Animal = hM4Di_animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage '\MovAnalysis\' Animal '_CuedMov_All.mat']);
    traces_files = dir(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace']);
    traces_files = {traces_files.name};
    traces_files = traces_files(3:end);
    for ii = 1:2
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep traces_files{ii} filesep...
            Animal '_' traces_files{ii} '_EphusTraces.mat'],'Lever_active','Lever_force_resample');
        temp_Bpod = [];
        responded_trials = [];
        temp_lever_active_session = [];
        clear speed_mvm_rwd;
        for jj = 1:length(CuedMov_SingleAnimal{ii}.Bpod)
            temp_Bpod = [temp_Bpod;CuedMov_SingleAnimal{ii}.Bpod{jj}];
            Rewarded_trial = find(~isnan(CuedMov_SingleAnimal{ii}.Bpod{jj}(:,5)));
            responded_trials = [responded_trials;unique([CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj}(:,1);Rewarded_trial])];
            session_end = nanmax(CuedMov_SingleAnimal{ii}.Bpod{jj}(end,:));
            session_end = round(session_end*1000);
            temp_lever_active_session = [temp_lever_active_session;Lever_active{jj}(1:session_end)];
            temp_info = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj};
            reward_index = ~isnan(temp_info(:,7));
            rewared_mvm_index = temp_info(:,5)==1;
            useful_index = logical(reward_index.*rewared_mvm_index);
            time_info = temp_info(useful_index,[2,7]);
            butterworth_stop = 5/500; % fraction of nyquist (cutoff = 10 Hz)
            [b a] = butter(4, butterworth_stop,'low');
            lever_force_smooth = filtfilt(b,a,Lever_force_resample{jj});
            lever_velocity_resample = [0;diff(lever_force_smooth)];
            lever_velocity_resample_smooth = smooth(lever_velocity_resample,1);
            time_info = round(time_info*1000);            
            for tt = 1:length(time_info)
                speed_mvm_rwd{jj,1}(tt,1) = nanmean(abs((lever_velocity_resample_smooth(time_info(tt,1):time_info(tt,2)))));
            end               
        end
        speed_mvm_rwd = cell2mat(speed_mvm_rwd);
        hM4Di_speed_mvm_rwd(curr_animal,ii) = nanmean(speed_mvm_rwd);
        trial_num = length(temp_Bpod);
        responded_trials = length(responded_trials);
        hM4Di_rspdfraction(curr_animal,ii) = responded_trials/trial_num;
        hM4Di_mvmfraction(curr_animal,ii) = sum(temp_lever_active_session)/length(temp_lever_active_session);
        hM4Di_rwd_in_all(curr_animal,ii) = size(CuedMov_SingleAnimal{ii}.LeverTrace_Reward_downsample,2)/size(CuedMov_SingleAnimal{ii}.LeverTrace_All_downsample,2);
        
    end
    clear CuedMov_SingleAnimal
end
        
for curr_animal = 1:length(mCherry_animals)
    Animal = mCherry_animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage '\MovAnalysis\'  Animal '_CuedMov_All.mat']);
    traces_files = dir(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace']);
    traces_files = {traces_files.name};
    traces_files = traces_files(3:end);
    for ii = 1:2
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep traces_files{ii} filesep...
            Animal '_' traces_files{ii} '_EphusTraces.mat'],'Lever_active','Lever_force_resample');
        temp_Bpod = [];
        responded_trials = [];
        temp_lever_active_session = [];
        clear speed_mvm_rwd;
        for jj = 1:length(CuedMov_SingleAnimal{ii}.Bpod)
            temp_Bpod = [temp_Bpod;CuedMov_SingleAnimal{ii}.Bpod{jj}];
            Rewarded_trial = find(~isnan(CuedMov_SingleAnimal{ii}.Bpod{jj}(:,5)));
            responded_trials = [responded_trials;unique([CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj}(:,1);Rewarded_trial])];
            session_end = nanmax(CuedMov_SingleAnimal{ii}.Bpod{jj}(end,:));
            session_end = round(session_end*1000);
            temp_lever_active_session = [temp_lever_active_session;Lever_active{jj}(1:session_end)];
            temp_info = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj};
            reward_index = ~isnan(temp_info(:,7));
            rewared_mvm_index = temp_info(:,5)==1;
            useful_index = logical(reward_index.*rewared_mvm_index);
            time_info = temp_info(useful_index,[2,7]);
            butterworth_stop = 5/500; % fraction of nyquist (cutoff = 10 Hz)
            [b a] = butter(4, butterworth_stop,'low');
            lever_force_smooth = filtfilt(b,a,Lever_force_resample{jj});
            lever_velocity_resample = [0;diff(lever_force_smooth)];
            lever_velocity_resample_smooth = smooth(lever_velocity_resample,1);
            time_info = round(time_info*1000);            
            for tt = 1:length(time_info)
                speed_mvm_rwd{jj,1}(tt,1) = nanmean(abs((lever_velocity_resample_smooth(time_info(tt,1):time_info(tt,2)))));
            end            
        end
        speed_mvm_rwd = cell2mat(speed_mvm_rwd);
        mCherry_speed_mvm_rwd(curr_animal,ii) = nanmean(speed_mvm_rwd);
        trial_num = length(temp_Bpod);
        responded_trials = length(responded_trials);
        mCherry_rspdfraction(curr_animal,ii) = responded_trials/trial_num;
        mCherry_mvmfraction(curr_animal,ii) = sum(temp_lever_active_session)/length(temp_lever_active_session);
        mCherry_rwd_in_all(curr_animal,ii) = size(CuedMov_SingleAnimal{ii}.LeverTrace_Reward_downsample,2)/size(CuedMov_SingleAnimal{ii}.LeverTrace_All_downsample,2);
        
    end
    clear CuedMov_SingleAnimal
end

save(['Z:\People\Chi\TwoP_IN\' IN '\' IN '_Behavior_Nai.mat'],'-v7.3');