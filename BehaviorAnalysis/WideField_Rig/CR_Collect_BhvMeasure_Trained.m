%% 
clearvars -except IN PV SOM VIP

% Initial = 'CR';
% Animals = {'3161016-O','3161018-R','3233232-L','3233232-O','3233232-R','3233233-O','3233233-L','3233233-R','3491479-L','3491479-LR','3547207-LR'};
% IN = 'PV';
% 
% Initial = 'CR';
% Animals = {'3183958-L','3183958-R','3183959-LL','3218181-L','3218181-O','3218181-R','3438521-O','3438521-L','3438521-R','3438521-LR','3453262-O','3453262-L'};
% IN = 'SOM';

Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};
IN = 'VIP';


for ii = 1:length(Animals)
    tic
    Animal = Animals{ii};
    disp(Animal)
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'ForBhv' filesep 'MovAnalysis']);
    load([Initial '_' Animal '_CuedMov.mat']);

    CR(ii,:) = nan(1,23);
    CRM_Dur(ii,:) = nan(1,23);
    Reaction_Time(ii,:) = nan(1,23);
    Cue_to_CuedRewardedMov(ii,:) = nan(1,23);
    CuedRewardedMov_to_Reward(ii,:) = nan(1,23);
    Cue_to_Reward(ii,:) = nan(1,23);
    MovSpeed(ii,:) = nan(1,23);
    TrialNum(ii,:) = nan(1,23);
    Corr_Reward_Matrix(:,:,ii)= nan(23,23);
    Corr_All_Matrix(:,:,ii)= nan(23,23);
    Corr_First_Matrix(:,:,ii)= nan(23,23);
    Corr_Reward_within(ii,:) = nan(1,23);
    Corr_Reward_across(ii,:) = nan(1,22);
    Corr_All_within(ii,:) = nan(1,23);
    Corr_All_across(ii,:) = nan(1,22);
    Corr_First_within(ii,:) = nan(1,23);
    Corr_First_across(ii,:) = nan(1,22);
    
    for jj = 1:length(CuedMov_SingleAnimal) %usefuldays
        if ~isempty(CuedMov_SingleAnimal{jj}) && ~isempty(CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All)
            CR(ii,jj) = CuedMov_SingleAnimal{jj}.CR;
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
            Cue_to_Reward(ii,jj) = nanmedian(temp(index_reward,7)-temp(index_reward,6));
            % Speed
            MovSpeed(ii,jj) = nanmean(nanmean(abs(diff(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample))));
            TrialNum(ii,jj) = size(CuedMov_SingleAnimal{jj}.LeverTrace_Reward_downsample,2);
        else
            CR(ii,jj) = nan;
            CRM_Dur(ii,jj) = nan;
            Reaction_Time(ii,jj) = nan;
            Cue_to_CuedRewardedMov(ii,jj) = nan;
            CuedRewardedMov_to_Reward(ii,jj) = nan;
            Cue_to_Reward(ii,jj) = nan;
            MovSpeed(ii,jj) = nan;
            TrialNum(ii,jj) = 0;
        end
        clear temp index_first index_reward
    end
        % Correlation
        Corr_Reward_Matrix(1:length(CuedMov_SingleAnimal),1:length(CuedMov_SingleAnimal),ii) = Trial_Trial_Corr_Reward;
        Corr_Reward_within(ii,:) = diag(Corr_Reward_Matrix(:,:,ii));
        Corr_Reward_across(ii,:) = diag(Corr_Reward_Matrix(:,:,ii),1);
        
        Corr_All_Matrix(1:length(CuedMov_SingleAnimal),1:length(CuedMov_SingleAnimal),ii) = Trial_Trial_Corr_All;
        Corr_All_within(ii,:) = diag(Corr_All_Matrix(:,:,ii));
        Corr_All_across(ii,:) = diag(Corr_All_Matrix(:,:,ii),1);
        
        Corr_First_Matrix(1:length(CuedMov_SingleAnimal),1:length(CuedMov_SingleAnimal),ii) = Trial_Trial_Corr_First;
        Corr_First_within(ii,:) = diag(Corr_First_Matrix(:,:,ii));
        Corr_First_across(ii,:) = diag(Corr_First_Matrix(:,:,ii),1);
        
        toc
end

MovSpeed = MovSpeed./10;

switch IN
    case 'PV'        
        PV.Animals = Animals;
        PV.CR = CR;
        PV.Reaction_Time = Reaction_Time;
        PV.MovSpeed = MovSpeed;
        PV.CRM_Dur = CRM_Dur;
        PV.Cue_to_CuedRewardedMov = Cue_to_CuedRewardedMov;
        PV.CuedRewardedMov_to_Reward = CuedRewardedMov_to_Reward;
        PV.Cue_to_Reward = Cue_to_Reward;
        PV.Corr_Reward_Matrix_mean = Corr_Reward_Matrix_mean;
        PV.Corr_All_Matrix_mean = Corr_All_Matrix_mean;
        PV.Corr_Reward_Matrix = Corr_Reward_Matrix;
        PV.Corr_All_Matrix = Corr_All_Matrix;
        PV.Corr_Reward_within = Corr_Reward_within;
        PV.Corr_Reward_across = Corr_Reward_across;
        PV.Corr_All_within = Corr_All_within;
        PV.Corr_All_across = Corr_All_across;
        PV.TrialNum = TrialNum;
    case 'SOM'        
        SOM.Animals = Animals;
        SOM.CR = CR;
        SOM.Reaction_Time = Reaction_Time;
        SOM.MovSpeed = MovSpeed;
        SOM.CRM_Dur = CRM_Dur;
        SOM.Cue_to_CuedRewardedMov = Cue_to_CuedRewardedMov;
        SOM.CuedRewardedMov_to_Reward = CuedRewardedMov_to_Reward;
        SOM.Cue_to_Reward = Cue_to_Reward;
        SOM.Corr_Reward_Matrix = Corr_Reward_Matrix;
        SOM.Corr_All_Matrix = Corr_All_Matrix;
        SOM.Corr_First_Matrix = Corr_First_Matrix;
        SOM.Corr_Reward_within = Corr_Reward_within;
        SOM.Corr_Reward_across = Corr_Reward_across;
        SOM.Corr_All_within = Corr_All_within;
        SOM.Corr_All_across = Corr_All_across;
        SOM.Corr_First_within = Corr_First_within;
        SOM.Corr_First_across = Corr_First_across;
        SOM.TrialNum = TrialNum;
    case 'VIP'        
        VIP.Animals = Animals;
        VIP.CR = CR;
        VIP.Reaction_Time = Reaction_Time;
        VIP.MovSpeed = MovSpeed;
        VIP.CRM_Dur = CRM_Dur;
        VIP.Cue_to_CuedRewardedMov = Cue_to_CuedRewardedMov;
        VIP.CuedRewardedMov_to_Reward = CuedRewardedMov_to_Reward;
        VIP.Cue_to_Reward = Cue_to_Reward;
        VIP.Corr_Reward_Matrix = Corr_Reward_Matrix;
        VIP.Corr_All_Matrix = Corr_All_Matrix;
        VIP.Corr_First_Matrix = Corr_First_Matrix;
        VIP.Corr_Reward_within = Corr_Reward_within;
        VIP.Corr_Reward_across = Corr_Reward_across;
        VIP.Corr_All_within = Corr_All_within;
        VIP.Corr_All_across = Corr_All_across;
        VIP.Corr_First_within = Corr_First_within;
        VIP.Corr_First_across = Corr_First_across;
        VIP.TrialNum = TrialNum;
end

cd(['Z:\People\Chi\WFLP_IN']);
save('All_BehaviorCollect_withPupilAnimal_refined_ForBhv.mat','VIP','SOM','-v7.3');
