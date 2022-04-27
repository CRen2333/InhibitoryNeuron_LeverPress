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

% Confirmed with immunostaining
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
    Corr_Reward_Matrix(:,:,ii) = nan(21,21);
    Trial_Num(ii,:) = nan(1,21);

    
    for jj = 1:length(CuedMov_SingleAnimal)%usefuldays
        if ~isempty(CuedMov_SingleAnimal{jj}) && ~isempty(CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All)
            % Cued rewarded movement Duration
            CRM_Dur(ii,jj) = CuedMov_SingleAnimal{jj}.Median_Movduration_Reward;
            CR(ii,jj) = CuedMov_SingleAnimal{jj}.CR;
            % Reaction time
            temp = CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All;
            index_first = temp(:,4)==1;
            Reaction_Time(ii,jj) = nanmedian(temp(index_first,2)-temp(index_first,6));
            % Cue to cued rewarded mov onset
            index_reward = temp(:,5)==1;
            Cue_to_CuedRewardedMov(ii,jj) = nanmedian(temp(index_reward,2)-temp(index_reward,6));
            CuedRewardedMov_to_Reward(ii,jj) = nanmedian(temp(index_reward,7)-temp(index_reward,2));
            Trial_Num(ii,jj) = size(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample,2); 
        else
            CRM_Dur(ii,jj) = nan;
            Reaction_Time(ii,jj) = nan;
            Cue_to_CuedRewardedMov(ii,jj) = nan;
            CuedRewardedMov_to_Reward(ii,jj) = nan;
            CR(ii,jj) = nan;
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

CR(3,1) = 0.19; % missing ephus file that day but still have dispatcher data

clear CuedMov_SingleAnimal
save('Z:\People\Chi\WFLP_IN\ChAT_Cas\Bhv_Collect.mat','-v7.3');


