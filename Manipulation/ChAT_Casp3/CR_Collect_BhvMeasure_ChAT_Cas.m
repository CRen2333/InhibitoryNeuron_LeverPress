%% 
close all
clear all
clc

Initial = 'CR';
IN = 'ChAT_Cas';
Animals = {'3702618-O','3702618-L','3702169-O','3702169-L','3702169-R','3702169-LR','3886988-O','3886988-L','3886989-L','3886990-O','3886990-L',...
    '3886990-R','3919112-O','3919112-L','3919112-R','3919114-O','3919114-L','3936512-O','3936512-L','3936513-O','3936513-L','3974586-LR','3974586-LL',...
    '3974585-O','3974585-L','3974585-R','3974584-L','4017416-O','4017416-L','4017417-O','4017417-L','4017417-R','4079574-L','4079632-L','4079632-O',...
    '4079633-O','4079633-L','4079633-R'};

Cas_Animals = {'3702618-O','3702169-O','3702169-L','3886988-O','3886989-L','3886990-O','3919112-O',...
    '3919112-L','3919114-O','3936512-O','3936513-O','3974586-LR','3974585-O','3974585-L','4017416-O','4017417-O','4079574-L','4079632-L','4079633-O'};
Sal_Animals = {'3702618-L','3702169-R','3702169-LR','3886988-L','3886990-L','3886990-R','3919112-R',...
    '3919114-L','3936512-L','3936513-L','3974586-LL','3974585-R','3974584-L','4017416-L','4017417-L','4017417-R','4079632-O','4079633-L','4079633-R'};

% Confirmed based on immunostaining
Cas_Animals_good = {'3702618-O','3702169-O','3702169-L','3886988-O','3886989-L','3919112-O','3919112-L',...
    '3919114-O','3936512-O','3936513-O','3974586-LR','3974585-O','3974585-L','4017416-O','4017417-O','4079574-L','4079632-L','4079633-O'};
Sal_Animals_good = {'3702618-L','3702169-R','3702169-LR','3886988-L','3886990-L','3886990-R','3919114-L',...
    '3936513-L','3974586-LL','3974584-L','4017416-L','4017417-L','4017417-R','4079632-O','4079633-L','4079633-R'};

Cas_index_good = ismember(Animals,Cas_Animals_good);
Sal_index_good = ismember(Animals,Sal_Animals_good);
Animals_inuse = [Cas_Animals_good,Sal_Animals_good];
inuse_index = ismember(Animals,Animals_inuse);

for ii = 1:length(Animals)
    Cages{ii} = Animals{ii}(1:7);
end

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
    MovSpeed(ii,:) = nan(1,21);
    Corr_Reward_Matrix(:,:,ii) = nan(21,21);
    Corr_All_Matrix(:,:,ii) = nan(21,21);
    PC1_fraction(ii,:) = nan(1,21);
    PC1_fraction_first(ii,:) = nan(1,21);
    PC1_fraction_all(ii,:) = nan(1,21);
    Reaction_Time_var(ii,:) = nan(1,21);
    Cue_to_CuedRewardedMov_var(ii,:) = nan(1,21);
    CuedRewardedMov_to_Reward_var(ii,:) = nan(1,21);
    Trial_Num(ii,:) = nan(1,21);
    Mov_order(ii,:) = nan(1,21);
    Mov_order_mean(ii,:) = nan(1,21);
    Mov_order_var(ii,:) = nan(1,21);
    RMS_reward(ii,:) = nan(1,21);
    RMS_all(ii,:) = nan(1,21);
    RMS_first(ii,:) = nan(1,21);
    
    for jj = 1:length(CuedMov_SingleAnimal)%usefuldays
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
            % Speed
            MovSpeed(ii,jj) = nanmean(nanmean(abs(diff(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample))))./10;
            [~,~,latent] = pca(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample);
            PC1_fraction(ii,jj) = latent(1)./sum(latent);
            [~,~,latent] = pca(CuedMov_SingleAnimal{jj}.LeverTrace_First_downsample);
            PC1_fraction_first(ii,jj) = latent(1)./sum(latent);
            [~,~,latent] = pca(CuedMov_SingleAnimal{jj}.LeverTrace_All_downsample);
            PC1_fraction_all(ii,jj) = latent(1)./sum(latent);
            Trial_Num(ii,jj) = size(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample,2);
            % RMS
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_reward(ii,jj) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_All_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_all(ii,jj) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_First_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_first(ii,jj) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
            temp_trace = CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample;
            temp_trace = temp_trace-repmat(temp_trace(1,:),201,1);
            MaxDisp(ii,jj) = nanmean(max(abs(temp_trace),[],1));
        else
            CRM_Dur(ii,jj) = nan;
            Reaction_Time(ii,jj) = nan;
            Cue_to_CuedRewardedMov(ii,jj) = nan;
            CuedRewardedMov_to_Reward(ii,jj) = nan;
            MovSpeed(ii,jj) = nan;
            PC1_fraction(ii,jj) = nan;
            CR(ii,jj) = nan;
            Reaction_Time_var(ii,jj) = nan;
            Cue_to_CuedRewardedMov_var(ii,jj) = nan;
            CuedRewardedMov_to_Reward_var(ii,jj) = nan;
            Mov_order(ii,jj) = nan;
            Mov_order_mean(ii,jj) = nan;
            Mov_order_var(ii,jj) = nan;
            RMS_reward(ii,jj) = nan;
            RMS_all(ii,jj) = nan;
            RMS_first(ii,jj) = nan;
            MaxDisp(ii,jj) = nan;
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
CR(3,1) = 0.19;
clear CuedMov_SingleAnimal
save('Z:\People\Chi\WFLP_IN\ChAT_Cas\Bhv_Collect.mat','-v7.3');


