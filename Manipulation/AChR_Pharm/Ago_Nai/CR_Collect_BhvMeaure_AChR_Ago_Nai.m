%% 
close all
clear all
clc

Initial = 'CR';
IN = 'AChR_Ago_Nai';
Animals = {'4259757-L','4302983-O','4259719-R','4302953-L','4333447-O','4333446-L','4383112-O','4383112-L','4383112-R','4383112-LR'};

Ago_Animals = {'4302983-O','4259719-R','4333447-O','4383112-R','4383112-LR'};
Sal_Animals = {'4259757-L','4302953-L','4333446-L','4383112-O','4383112-L'};

Ago_Animals_good = {'4302983-O','4259719-R','4333447-O','4383112-R','4383112-LR'};
Sal_Animals_good = {'4259757-L','4302953-L','4333446-L','4383112-O','4383112-L'};

Ago_index = ismember(Animals,Ago_Animals);
Sal_index = ismember(Animals,Sal_Animals);
Ago_index_good = ismember(Animals,Ago_Animals_good);
Sal_index_good = ismember(Animals,Sal_Animals_good);
Animals_inuse = [Ago_Animals_good,Sal_Animals_good];
inuse_index = ismember(Animals,Animals_inuse);

for ii = 1:length(Animals)
    tic
    Animal = Animals{ii};
    disp(Animal)
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'MovAnalysis']);
    load([Initial '_' Animal '_CuedMov.mat']);

    usefuldays = [1:4];
    CR(ii,:) = nan(1,4);
    CRM_Dur(ii,:) = nan(1,4);
    Reaction_Time(ii,:) = nan(1,4);
    Cue_to_CuedRewardedMov(ii,:) = nan(1,4);
    CuedRewardedMov_to_Reward(ii,:) = nan(1,4);
    MovSpeed(ii,:) = nan(1,4);
    Corr_Reward_Matrix(:,:,ii) = nan(4,4);
    Corr_All_Matrix(:,:,ii) = nan(4,4);
    Reaction_Time_var(ii,:) = nan(1,4);
    Cue_to_CuedRewardedMov_var(ii,:) = nan(1,4);
    CuedRewardedMov_to_Reward_var(ii,:) = nan(1,4);
    Trial_Num(ii,:) = nan(1,4);
    MaxDisp(ii,:) = nan(1,4);
        
    for jj = 1:length(CuedMov_SingleAnimal) %usefuldays
        if ~isempty(CuedMov_SingleAnimal{jj}) && ~isempty(CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All)
            % Cued rewarded movement Duration
            CRM_Dur(ii,jj) = CuedMov_SingleAnimal{jj}.Median_Movduration_Reward;
            CR(ii,jj) = CuedMov_SingleAnimal{jj}.CR;
            % Reaction time
            temp = CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All;
            index_first = temp(:,4)==1;
            Reaction_Time(ii,jj) = nanmedian(temp(index_first,2)-temp(index_first,6));
            Reaction_Time_var(ii,jj) = std(temp(index_first,2)-temp(index_first,6));
            % Cue to cued rewarded mov onset
            index_reward = temp(:,5)==1;
            Mov_order(ii,jj) = nanmedian(temp(index_reward,4));
            Mov_order_mean(ii,jj) = nanmean(temp(index_reward,4));
            Mov_order_var(ii,jj) = nanstd(temp(index_reward,4));
            Cue_to_CuedRewardedMov(ii,jj) = nanmedian(temp(index_reward,2)-temp(index_reward,6));
            CuedRewardedMov_to_Reward(ii,jj) = nanmedian(temp(index_reward,7)-temp(index_reward,2));
            Cue_to_CuedRewardedMov_var(ii,jj) = std(temp(index_reward,2)-temp(index_reward,6));
            CuedRewardedMov_to_Reward_var(ii,jj) = std(temp(index_reward,7)-temp(index_reward,2));
            Trial_Num(ii,jj) = size(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample,2);

            % Speed & Max displacement
            MovSpeed(ii,jj) = nanmean(nanmean(abs(diff(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample))))./10;
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample;
            temp_trace = temp_trace-repmat(temp_trace(1,:),201,1);
            MaxDisp(ii,jj) = nanmean(max(abs(temp_trace),[],1));
            
            % fraction of give up trials
            
        end               
        clear temp index_first index_reward temp_trace
    end
        % Correlation
        days = size(Trial_Trial_Corr_Reward,1);
        Corr_Reward_Matrix(1:days,1:days,ii) = Trial_Trial_Corr_Reward;
        Corr_Reward_within(ii,:) = diag(Corr_Reward_Matrix(:,:,ii));
        Corr_Reward_across(ii,:) = diag(Corr_Reward_Matrix(:,:,ii),1);
        
        Corr_All_Matrix(1:days,1:days,ii) = Trial_Trial_Corr_All;
        Corr_All_within(ii,:) = diag(Corr_All_Matrix(:,:,ii));
        Corr_All_across(ii,:) = diag(Corr_All_Matrix(:,:,ii),1);        
        
        Corr_First_Matrix(1:days,1:days,ii) = Trial_Trial_Corr_First;
        Corr_First_within(ii,:) = diag(Corr_First_Matrix(:,:,ii));
        Corr_First_across(ii,:) = diag(Corr_First_Matrix(:,:,ii),1);        
        toc
end
clear CuedMov_SingleAnimal

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\WFLP_IN\CheckNonconsTrials' filesep IN '_' Initial '_' Animal '_CheckNonconsTrials_Corrected.mat'], 'TrialInfo_xsg_all');
    load(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\MovAnalysis\' Initial '_' Animal '_Traces.mat'],'Lever_Active');
    for ii = 1:4
        temp_lever_active_session = [];
        TrialInfo_xsg_all{ii} = TrialInfo_xsg_all{ii}(~cellfun(@isempty, TrialInfo_xsg_all{ii}));
        for jj = 1:min(length(TrialInfo_xsg_all{ii}),length(Lever_Active{ii}))
            if jj == min(length(TrialInfo_xsg_all{ii}),length(Lever_Active{ii}))
                session_end = TrialInfo_xsg_all{ii}{jj}(end,1);
            else
                session_end = 300;
            end
            session_end = round(session_end*1000);
            temp_lever_active_session = [temp_lever_active_session;Lever_Active{ii}{jj}(1:session_end)];
        end
        mvmfraction(curr_animal,ii) = sum(temp_lever_active_session)/length(temp_lever_active_session);

    end
    clear Lever_Active TrialInfo_xsg_all
end

save(['Z:\People\Chi\WFLP_IN\' IN '\Bhv_Collect.mat'],'-v7.3');

