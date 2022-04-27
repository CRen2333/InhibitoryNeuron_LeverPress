%% 
close all
clear all
clc

Initial = 'CR';
IN = 'ChAT_eOPN3';
Animals = {'4333446-O','4333446-R','4333446-LR','4365811-LL','4365811-O','4365811-L','4365811-R','4383177-L','4383177-R','4383177-LR','4412583-O','4412583-L','4412583-R','4412582-O','4412582-L','4412601-O','4412601-L'};

eOPN3_Animals = {'4333446-O','4333446-R','4365811-O','4383177-L','4383177-R','4412583-O','4412583-L','4412583-R','4412601-L'};
mCherry_Animals = {'4333446-LR','CR_4365811-LL','4365811-L','4365811-R','4383177-LR','4412582-O','4412582-L','4412601-O'};

% Confirmed with histology
eOPN3_Animals_good = {'4333446-O','4333446-R','4365811-O','4383177-L','4383177-R','4412583-O','4412583-L','4412583-R','4412601-L'};
mCherry_Animals_good = {'4333446-LR','4365811-LL','4365811-L','4365811-R','4383177-LR','4412582-O','4412582-L','4412601-O'};

eOPN3_index = ismember(Animals,eOPN3_Animals);
mCherry_index = ismember(Animals,mCherry_Animals);
eOPN3_index_good = ismember(Animals,eOPN3_Animals_good);
mCherry_index_good = ismember(Animals,mCherry_Animals_good);
Animals_inuse = [eOPN3_Animals_good,mCherry_Animals_good];
inuse_index = ismember(Animals,Animals_inuse);

for ii = 1:length(Animals)
    tic
    Animal = Animals{ii};
    disp(Animal)
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'MovAnalysis']);
    load([Initial '_' Animal '_CuedMov.mat']);

    usefuldays = [1:21];
    CR(ii,:) = nan(1,21);
    CRM_Dur(ii,:) = nan(1,21);
    Reaction_Time(ii,:) = nan(1,21);
    Cue_to_CuedRewardedMov(ii,:) = nan(1,21);
    CuedRewardedMov_to_Reward(ii,:) = nan(1,21);
    Corr_Reward_Matrix(:,:,ii) = nan(21,21);
    Trial_Num(ii,:) = nan(1,21);
    MaxDisp(ii,:) = nan(1,21);
        
    for jj = 1:length(CuedMov_SingleAnimal)%usefuldays
        if ~isempty(CuedMov_SingleAnimal{jj})
            CR(ii,jj) = CuedMov_SingleAnimal{jj}.CR;
        end
        if ~isempty(CuedMov_SingleAnimal{jj}) && ~isempty(CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All)
            % Cued rewarded movement Duration
            CRM_Dur(ii,jj) = CuedMov_SingleAnimal{jj}.Median_Movduration_Reward;
            % Reaction time
            temp = CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All;
            index_first = temp(:,4)==1;
            Reaction_Time(ii,jj) = nanmedian(temp(index_first,2)-temp(index_first,6));
            % Cue to cued rewarded mov onset
            index_reward = temp(:,5)==1;
            Cue_to_CuedRewardedMov(ii,jj) = nanmedian(temp(index_reward,2)-temp(index_reward,6));
            CuedRewardedMov_to_Reward(ii,jj) = nanmedian(temp(index_reward,7)-temp(index_reward,2));
            Trial_Num(ii,jj) = size(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample,2);

            % Max displacement
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample;
            temp_trace = temp_trace-repmat(temp_trace(1,:),201,1);
            MaxDisp(ii,jj) = nanmean(max(abs(temp_trace),[],1));           
          
        end               
        clear temp index_first index_reward temp_trace
    end
    % Correlation
    days = size(Trial_Trial_Corr_Reward,1);
    Corr_Reward_Matrix(1:days,1:days,ii) = Trial_Trial_Corr_Reward;
    Corr_Reward_within(ii,:) = diag(Corr_Reward_Matrix(:,:,ii));
    Corr_Reward_across(ii,:) = diag(Corr_Reward_Matrix(:,:,ii),1);

    toc
end
clear CuedMov_SingleAnimal

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\WFLP_IN\CheckNonconsTrials' filesep IN '_' Initial '_' Animal '_CheckNonconsTrials_Corrected.mat'], 'TrialInfo_xsg_all');
    load(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\MovAnalysis\' Initial '_' Animal '_Traces.mat'],'Lever_Active');
    for ii = 1:length(TrialInfo_xsg_all)
        if isempty(TrialInfo_xsg_all{ii})
            mvmfraction(curr_animal,ii) = nan;
            continue
        end

        temp_lever_active_session = [];
        temp_index = ~cellfun(@isempty, TrialInfo_xsg_all{ii});
        temp_index = temp_index(1:min(length(TrialInfo_xsg_all{ii}),length(Lever_Active{ii})));
        TrialInfo_xsg_all{ii} = TrialInfo_xsg_all{ii}(temp_index);
        Lever_Active{ii} = Lever_Active{ii}(temp_index);
        
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
cd(['Z:\People\Chi\WFLP_IN\' IN]);

