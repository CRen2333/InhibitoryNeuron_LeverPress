%% Collect_behavior_2P
clearvars -except SOM_2P VIP_2P ChAT_2P;
close all;
clc;

IN = 'SOM';
switch IN
    case 'VIP'
        Animals = {'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','WL_3526641-R','WL_3526642-L','WL_3526642-R'};
    case 'SOM'
        Animals = {'KP_3459921_1','KP_3461990_1','WL_3526578-O','WL_3526580-O','WL_3547273-LR','WL_3547273-R','CR_3786142-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'};
    case 'ChAT'
        Animals = {'CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R'};
    case 'VIP_SOM'
        Animals = {'CR_3974574-R','CR_3974573-O','CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-L','CR_4017421-R',...
            'CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042830-R','CR_4042832-O','CR_4042832-L','CR_4042832-R'};
    case 'VIP_hM4Di'
        Animals = {'CR_4113793-O','CR_4113793-L'};
end

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal '\MovAnalysis\' Animal '_CuedMov_All.mat']);
    temp_matrix = nan(22,22);
    temp_matrix(1:size(Trial_Trial_Corr_Reward,1),1:size(Trial_Trial_Corr_Reward,1)) = Trial_Trial_Corr_Reward;
    LeverCorr_Reward_Matrix(:,:,curr_animal) = temp_matrix;
    temp_matrix = nan(22,22);
    temp_matrix(1:size(Trial_Trial_Corr_All,1),1:size(Trial_Trial_Corr_All,1)) = Trial_Trial_Corr_All;
    LeverCorr_All_Matrix(:,:,curr_animal) = temp_matrix;
    temp_matrix = nan(22,22);
    temp_matrix(1:size(Trial_Trial_Corr_First,1),1:size(Trial_Trial_Corr_First,1)) = Trial_Trial_Corr_First;
    LeverCorr_First_Matrix(:,:,curr_animal) = temp_matrix;
    RwdMov_Duration(curr_animal,[1:22]) = nan;
    CR(curr_animal,[1:22]) = nan;
    RT(curr_animal,[1:22]) = nan;
    RT_var(curr_animal,[1:22]) = nan;
    RwdMVM_Rwd(curr_animal,[1:22]) = nan;
    RwdMVM_order(curr_animal,[1:22]) = nan;
    RwdMVM_order_mean(curr_animal,[1:22]) = nan;
    C_CRM(curr_animal,[1:22]) = nan;
    C_Rwd(curr_animal,[1:22]) = nan;
    First_R(curr_animal,[1:22]) = nan;
    TrialNum(curr_animal,[1:22]) = nan;
    RMS_reward(curr_animal,[1:22]) = nan;
    RMS_all(curr_animal,[1:22]) = nan;
    RMS_first(curr_animal,[1:22]) = nan;
    RwdMVM_speed(curr_animal,[1:22]) = nan;
    FirstMVM_speed(curr_animal,[1:22]) = nan;
    AllMVM_speed(curr_animal,[1:22]) = nan;
    
    for curr_day = 1:length(CuedMov_SingleAnimal)
        if isempty(CuedMov_SingleAnimal{curr_day})
            continue
        end
        RwdMov_Duration(curr_animal,curr_day) = CuedMov_SingleAnimal{curr_day}.Median_Movduration_Reward;
        if strcmp(IN,'ChAT')
            CR(curr_animal,curr_day) = CuedMov_SingleAnimal{curr_day}.CR;
        else
            temp_bpod = cell2mat(CuedMov_SingleAnimal{curr_day}.Bpod');
            CR(curr_animal,curr_day) = sum(~isnan(temp_bpod(:,5)))./size(temp_bpod,1);
        end
        temp_info = CuedMov_SingleAnimal{curr_day}.Cued_MovOnset_Info_All;
        first_index = temp_info(:,4)==1;
        reward_index = ~isnan(temp_info(:,7));
        rewared_mvm_index = temp_info(:,5)==1;
        RT(curr_animal,curr_day) = nanmedian(temp_info(first_index,2)-temp_info(first_index,6));
        RT_var(curr_animal,curr_day) = nanstd(temp_info(first_index,2)-temp_info(first_index,6));
        RwdMVM_Rwd(curr_animal,curr_day) = nanmedian(temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2));
        C_CRM(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,2)-temp_info(rewared_mvm_index,6));
        C_Rwd(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,7)-temp_info(rewared_mvm_index,6));
        First_R(curr_animal,curr_day) = nanmedian(temp_info(first_index,7)-temp_info(first_index,2));
        RwdMVM_order(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,4));
        RwdMVM_order_mean(curr_animal,curr_day) = nanmean(temp_info(rewared_mvm_index,4));
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
        RwdMVM_speed(curr_animal,curr_day) = nanmean(nanmean(abs(diff(temp_trace,[],1)),1));
        TrialNum(curr_animal,curr_day) = size(temp_trace,2);
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_All_downsample;
        AllMVM_speed = nanmean(nanmean(abs(diff(temp_trace,[],1)),1));
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_First_downsample;
        FirstMVM_speed = nanmean(nanmean(abs(diff(temp_trace,[],1)),1));
        clear temp_trace
        % RMS
        if ~isempty(CuedMov_SingleAnimal{curr_day}) && ~isempty(CuedMov_SingleAnimal{curr_day}.Cued_MovOnset_Info_All)
            temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_reward(curr_animal,curr_day) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
            temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_All_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_all(curr_animal,curr_day) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
            temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_First_downsample;
            temp_trace_mean = nanmean(temp_trace,2);
            temp_trace = temp_trace-repmat(temp_trace_mean,1,size(temp_trace,2));
            temp_trace = temp_trace.*temp_trace;
            temp_trace = sqrt(sum(temp_trace,1)./size(temp_trace,1));
            RMS_first(curr_animal,curr_day) = sqrt(sum(temp_trace))/length(temp_trace)/sum(abs(diff(temp_trace_mean)));
        else
            RMS_reward(curr_animal,curr_day) = nan;
            RMS_all(curr_animal,curr_day) = nan;
            RMS_first(curr_animal,curr_day) = nan;
        end
    end
    clear clear temp_trace Trial_Trial_Corr_Reward Trial_Trial_Corr_First Trial_Trial_Corr_All CuedMov_SingleAnimal
end
switch IN
    case 'SOM'
        SOM_2P.LeverCorr_Reward_Matrix = LeverCorr_Reward_Matrix;
        SOM_2P.LeverCorr_All_Matrix = LeverCorr_All_Matrix;
        SOM_2P.LeverCorr_First_Matrix = LeverCorr_First_Matrix;
        SOM_2P.RwdMov_Duration = RwdMov_Duration;
        SOM_2P.CR = CR;
        SOM_2P.RT = RT;
        SOM_2P.RT_var = RT_var;
        SOM_2P.RwdMVM_Rwd = RwdMVM_Rwd;
        SOM_2P.RwdMVM_speed = RwdMVM_speed/10;
        SOM_2P.C_CRM = C_CRM;
        SOM_2P.C_Rwd = C_Rwd;
        SOM_2P.TrialNum = TrialNum;
    case 'VIP'
        VIP_2P.LeverCorr_Reward_Matrix = LeverCorr_Reward_Matrix;
        VIP_2P.LeverCorr_All_Matrix = LeverCorr_All_Matrix;
        VIP_2P.LeverCorr_First_Matrix = LeverCorr_First_Matrix;
        VIP_2P.RwdMov_Duration = RwdMov_Duration;
        VIP_2P.CR = CR;
        VIP_2P.RT = RT;
        VIP_2P.RT_var = RT_var;
        VIP_2P.RwdMVM_Rwd = RwdMVM_Rwd;
        VIP_2P.RwdMVM_speed = RwdMVM_speed/10;
        VIP_2P.C_CRM = C_CRM;
        VIP_2P.C_Rwd = C_Rwd;
        VIP_2P.TrialNum = TrialNum;
    case 'ChAT'
        ChAT_2P.LeverCorr_Reward_Matrix = LeverCorr_Reward_Matrix;
        ChAT_2P.LeverCorr_All_Matrix = LeverCorr_All_Matrix;
        ChAT_2P.LeverCorr_First_Matrix = LeverCorr_First_Matrix;
        ChAT_2P.RwdMov_Duration = RwdMov_Duration;
        ChAT_2P.CR = CR;
        ChAT_2P.RT = RT;
        ChAT_2P.RT_var = RT_var;
        ChAT_2P.RwdMVM_Rwd = RwdMVM_Rwd;
        ChAT_2P.PC1_Fraction = PC1_fraction;
        ChAT_2P.PC1_Score = PC1_Score;
        ChAT_2P.TrialNum = TrialNum;
    case 'VIP_SOM'
        VIPSOM_2P.LeverCorr_Reward_Matrix = LeverCorr_Reward_Matrix;
        VIPSOM_2P.LeverCorr_All_Matrix = LeverCorr_All_Matrix;
        VIPSOM_2P.LeverCorr_First_Matrix = LeverCorr_First_Matrix;
        VIPSOM_2P.RwdMov_Duration = RwdMov_Duration;
        VIPSOM_2P.CR = CR;
        VIPSOM_2P.RT = RT;
        VIPSOM_2P.RT_var = RT_var;
        VIPSOM_2P.RwdMVM_Rwd = RwdMVM_Rwd;
        VIPSOM_2P.C_CRM = C_CRM;
        VIPSOM_2P.First_R = First_R;
        VIPSOM_2P.PC1_Fraction = PC1_fraction;
        VIPSOM_2P.PC1_Fraction_time = PC1_fraction_time;
        VIPSOM_2P.PC1_Fraction_trial = PC1_fraction_trial;
        VIPSOM_2P.TrialNum = TrialNum;
        VIPSOM_2P.RMS_reward = RMS_reward;
        VIPSOM_2P.RMS_all = RMS_all;
        VIPSOM_2P.RMS_first = RMS_first;
end

%% 4 stages
trialnum_thre = 3;
VIP_2P.CR_stages(:,1) = nanmean(VIP_2P.CR(:,1:2),2);
VIP_2P.CR_stages(:,2) = nanmean(VIP_2P.CR(:,3:8),2);
VIP_2P.CR_stages(:,3) = nanmean(VIP_2P.CR(:,9:16),2);
VIP_2P.CR_stages(:,4) = nanmean(VIP_2P.CR(:,17:22),2);
temp_var = VIP_2P.RwdMov_Duration;
temp_var(VIP_2P.TrialNum<trialnum_thre) = nan;
VIP_2P.RwdMov_Duration_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.RwdMov_Duration_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.RwdMov_Duration_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.RwdMov_Duration_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.RT;
VIP_2P.RT_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.RT_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.RT_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.RT_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.RT_var;
VIP_2P.RT_var_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.RT_var_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.RT_var_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.RT_var_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.RwdMVM_Rwd;
temp_var(VIP_2P.TrialNum<trialnum_thre) = nan;
VIP_2P.RwdMVM_Rwd_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.RwdMVM_Rwd_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.RwdMVM_Rwd_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.RwdMVM_Rwd_stages(:,4) = nanmean(temp_var(:,17:22),2);
for ii = 1:size(VIP_2P.LeverCorr_Reward_Matrix,3)
    VIP_2P.LeverCorr_Reward_within(ii,:) = diag(VIP_2P.LeverCorr_Reward_Matrix(:,:,ii),0);
    VIP_2P.LeverCorr_Reward_across(ii,:) = diag(VIP_2P.LeverCorr_Reward_Matrix(:,:,ii),1);
    VIP_2P.LeverCorr_First_within(ii,:) = diag(VIP_2P.LeverCorr_First_Matrix(:,:,ii),0);
    VIP_2P.LeverCorr_First_across(ii,:) = diag(VIP_2P.LeverCorr_First_Matrix(:,:,ii),1);
    VIP_2P.LeverCorr_All_within(ii,:) = diag(VIP_2P.LeverCorr_All_Matrix(:,:,ii),0);
    VIP_2P.LeverCorr_All_across(ii,:) = diag(VIP_2P.LeverCorr_All_Matrix(:,:,ii),1);
end
temp_var = VIP_2P.LeverCorr_Reward_within;
temp_var(VIP_2P.TrialNum<trialnum_thre) = nan;
VIP_2P.LeverCorr_within_Reward_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.LeverCorr_within_Reward_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.LeverCorr_within_Reward_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.LeverCorr_within_Reward_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.LeverCorr_Reward_across;
trialnum_index_2 = (VIP_2P.TrialNum(:,1:end-1)+VIP_2P.TrialNum(:,2:end))<trialnum_thre;
temp_var(trialnum_index_2) = nan;
VIP_2P.LeverCorr_across_Reward_stages(:,1) = nanmean(temp_var(:,1),2);
VIP_2P.LeverCorr_across_Reward_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIP_2P.LeverCorr_across_Reward_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIP_2P.LeverCorr_across_Reward_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = VIP_2P.LeverCorr_First_within;
temp_var(VIP_2P.TrialNum<trialnum_thre) = nan;
VIP_2P.LeverCorr_within_First_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.LeverCorr_within_First_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.LeverCorr_within_First_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.LeverCorr_within_First_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.LeverCorr_First_across;
temp_var(trialnum_index_2) = nan;
VIP_2P.LeverCorr_across_First_stages(:,1) = nanmean(temp_var(:,1),2);
VIP_2P.LeverCorr_across_First_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIP_2P.LeverCorr_across_First_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIP_2P.LeverCorr_across_First_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = VIP_2P.LeverCorr_All_within;
temp_var(VIP_2P.TrialNum<trialnum_thre) = nan;
VIP_2P.LeverCorr_within_All_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIP_2P.LeverCorr_within_All_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIP_2P.LeverCorr_within_All_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIP_2P.LeverCorr_within_All_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIP_2P.LeverCorr_All_across;
temp_var(trialnum_index_2) = nan;
VIP_2P.LeverCorr_across_All_stages(:,1) = nanmean(temp_var(:,1),2);
VIP_2P.LeverCorr_across_All_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIP_2P.LeverCorr_across_All_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIP_2P.LeverCorr_across_All_stages(:,4) = nanmean(temp_var(:,17:21),2);

ChAT_2P.CR_stages(:,1) = nanmean(ChAT_2P.CR(:,1:2),2);
ChAT_2P.CR_stages(:,2) = nanmean(ChAT_2P.CR(:,3:8),2);
ChAT_2P.CR_stages(:,3) = nanmean(ChAT_2P.CR(:,9:16),2);
ChAT_2P.CR_stages(:,4) = nanmean(ChAT_2P.CR(:,17:22),2);
temp_var = ChAT_2P.RwdMov_Duration;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.RwdMov_Duration_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.RwdMov_Duration_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.RwdMov_Duration_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.RwdMov_Duration_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.RT;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.RT_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.RT_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.RT_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.RT_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.RT_var;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.RT_var_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.RT_var_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.RT_var_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.RT_var_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.RwdMVM_Rwd;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.RwdMVM_Rwd_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.RwdMVM_Rwd_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.RwdMVM_Rwd_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.RwdMVM_Rwd_stages(:,4) = nanmean(temp_var(:,17:22),2);
for ii = 1:size(ChAT_2P.LeverCorr_Reward_Matrix,3)
    ChAT_2P.LeverCorr_Reward_within(ii,:) = diag(ChAT_2P.LeverCorr_Reward_Matrix(:,:,ii),0);
    ChAT_2P.LeverCorr_Reward_across(ii,:) = diag(ChAT_2P.LeverCorr_Reward_Matrix(:,:,ii),1);
    ChAT_2P.LeverCorr_First_within(ii,:) = diag(ChAT_2P.LeverCorr_First_Matrix(:,:,ii),0);
    ChAT_2P.LeverCorr_First_across(ii,:) = diag(ChAT_2P.LeverCorr_First_Matrix(:,:,ii),1);
    ChAT_2P.LeverCorr_All_within(ii,:) = diag(ChAT_2P.LeverCorr_All_Matrix(:,:,ii),0);
    ChAT_2P.LeverCorr_All_across(ii,:) = diag(ChAT_2P.LeverCorr_All_Matrix(:,:,ii),1);
end
temp_var = ChAT_2P.LeverCorr_Reward_within;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.LeverCorr_within_Reward_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.LeverCorr_within_Reward_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.LeverCorr_within_Reward_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.LeverCorr_within_Reward_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.LeverCorr_Reward_across;
ChAT_2P.LeverCorr_across_Reward_stages(:,1) = nanmean(temp_var(:,1),2);
ChAT_2P.LeverCorr_across_Reward_stages(:,2) = nanmean(temp_var(:,3:7),2);
ChAT_2P.LeverCorr_across_Reward_stages(:,3) = nanmean(temp_var(:,9:15),2);
ChAT_2P.LeverCorr_across_Reward_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = ChAT_2P.LeverCorr_First_within;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.LeverCorr_within_First_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.LeverCorr_within_First_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.LeverCorr_within_First_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.LeverCorr_within_First_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.LeverCorr_First_across;
ChAT_2P.LeverCorr_across_First_stages(:,1) = nanmean(temp_var(:,1),2);
ChAT_2P.LeverCorr_across_First_stages(:,2) = nanmean(temp_var(:,3:7),2);
ChAT_2P.LeverCorr_across_First_stages(:,3) = nanmean(temp_var(:,9:15),2);
ChAT_2P.LeverCorr_across_First_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = ChAT_2P.LeverCorr_All_within;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.LeverCorr_within_All_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.LeverCorr_within_All_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.LeverCorr_within_All_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.LeverCorr_within_All_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = ChAT_2P.LeverCorr_All_across;
ChAT_2P.LeverCorr_across_All_stages(:,1) = nanmean(temp_var(:,1),2);
ChAT_2P.LeverCorr_across_All_stages(:,2) = nanmean(temp_var(:,3:7),2);
ChAT_2P.LeverCorr_across_All_stages(:,3) = nanmean(temp_var(:,9:15),2);
ChAT_2P.LeverCorr_across_All_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = ChAT_2P.PC1_Fraction;
temp_var(ChAT_2P.TrialNum<trialnum_thre) = nan;
ChAT_2P.PC1_Fraction_stages(:,1) = nanmean(temp_var(:,1:2),2);
ChAT_2P.PC1_Fraction_stages(:,2) = nanmean(temp_var(:,3:8),2);
ChAT_2P.PC1_Fraction_stages(:,3) = nanmean(temp_var(:,9:16),2);
ChAT_2P.PC1_Fraction_stages(:,4) = nanmean(temp_var(:,17:22),2);

SOM_2P.CR_stages(:,1) = nanmean(SOM_2P.CR(:,1:2),2);
SOM_2P.CR_stages(:,2) = nanmean(SOM_2P.CR(:,3:8),2);
SOM_2P.CR_stages(:,3) = nanmean(SOM_2P.CR(:,9:16),2);
SOM_2P.CR_stages(:,4) = nanmean(SOM_2P.CR(:,17:22),2);
temp_var = SOM_2P.RwdMov_Duration;
temp_var(SOM_2P.TrialNum<trialnum_thre) = nan;
SOM_2P.RwdMov_Duration_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.RwdMov_Duration_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.RwdMov_Duration_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.RwdMov_Duration_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.RT;
SOM_2P.RT_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.RT_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.RT_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.RT_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.RT_var;
SOM_2P.RT_var_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.RT_var_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.RT_var_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.RT_var_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.RwdMVM_Rwd;
temp_var(SOM_2P.TrialNum<trialnum_thre) = nan;
SOM_2P.RwdMVM_Rwd_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.RwdMVM_Rwd_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.RwdMVM_Rwd_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.RwdMVM_Rwd_stages(:,4) = nanmean(temp_var(:,17:22),2);
for ii = 1:size(SOM_2P.LeverCorr_Reward_Matrix,3)
    SOM_2P.LeverCorr_Reward_within(ii,:) = diag(SOM_2P.LeverCorr_Reward_Matrix(:,:,ii),0);
    SOM_2P.LeverCorr_Reward_across(ii,:) = diag(SOM_2P.LeverCorr_Reward_Matrix(:,:,ii),1);
    SOM_2P.LeverCorr_First_within(ii,:) = diag(SOM_2P.LeverCorr_First_Matrix(:,:,ii),0);
    SOM_2P.LeverCorr_First_across(ii,:) = diag(SOM_2P.LeverCorr_First_Matrix(:,:,ii),1);
    SOM_2P.LeverCorr_All_within(ii,:) = diag(SOM_2P.LeverCorr_All_Matrix(:,:,ii),0);
    SOM_2P.LeverCorr_All_across(ii,:) = diag(SOM_2P.LeverCorr_All_Matrix(:,:,ii),1);
end
temp_var = SOM_2P.LeverCorr_Reward_within;
temp_var(SOM_2P.TrialNum<trialnum_thre) = nan;
SOM_2P.LeverCorr_within_Reward_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.LeverCorr_within_Reward_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.LeverCorr_within_Reward_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.LeverCorr_within_Reward_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.LeverCorr_Reward_across;
trialnum_index_2 = (SOM_2P.TrialNum(:,1:end-1)+SOM_2P.TrialNum(:,2:end))<trialnum_thre;
temp_var(trialnum_index_2) = nan;
SOM_2P.LeverCorr_across_Reward_stages(:,1) = nanmean(temp_var(:,1),2);
SOM_2P.LeverCorr_across_Reward_stages(:,2) = nanmean(temp_var(:,3:7),2);
SOM_2P.LeverCorr_across_Reward_stages(:,3) = nanmean(temp_var(:,9:15),2);
SOM_2P.LeverCorr_across_Reward_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = SOM_2P.LeverCorr_First_within;
temp_var(SOM_2P.TrialNum<trialnum_thre) = nan;
SOM_2P.LeverCorr_within_First_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.LeverCorr_within_First_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.LeverCorr_within_First_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.LeverCorr_within_First_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.LeverCorr_First_across;
trialnum_index_2 = (SOM_2P.TrialNum(:,1:end-1)+SOM_2P.TrialNum(:,2:end))<trialnum_thre;
temp_var(trialnum_index_2) = nan;
SOM_2P.LeverCorr_across_First_stages(:,1) = nanmean(temp_var(:,1),2);
SOM_2P.LeverCorr_across_First_stages(:,2) = nanmean(temp_var(:,3:7),2);
SOM_2P.LeverCorr_across_First_stages(:,3) = nanmean(temp_var(:,9:15),2);
SOM_2P.LeverCorr_across_First_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = SOM_2P.LeverCorr_All_within;
temp_var(SOM_2P.TrialNum<trialnum_thre) = nan;
SOM_2P.LeverCorr_within_All_stages(:,1) = nanmean(temp_var(:,1:2),2);
SOM_2P.LeverCorr_within_All_stages(:,2) = nanmean(temp_var(:,3:8),2);
SOM_2P.LeverCorr_within_All_stages(:,3) = nanmean(temp_var(:,9:16),2);
SOM_2P.LeverCorr_within_All_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = SOM_2P.LeverCorr_All_across;
trialnum_index_2 = (SOM_2P.TrialNum(:,1:end-1)+SOM_2P.TrialNum(:,2:end))<trialnum_thre;
temp_var(trialnum_index_2) = nan;
SOM_2P.LeverCorr_across_All_stages(:,1) = nanmean(temp_var(:,1),2);
SOM_2P.LeverCorr_across_All_stages(:,2) = nanmean(temp_var(:,3:7),2);
SOM_2P.LeverCorr_across_All_stages(:,3) = nanmean(temp_var(:,9:15),2);
SOM_2P.LeverCorr_across_All_stages(:,4) = nanmean(temp_var(:,17:21),2);

trialnum_thre = 3;
VIPSOM_2P.CR_stages(:,1) = nanmean(VIPSOM_2P.CR(:,1:2),2);
VIPSOM_2P.CR_stages(:,2) = nanmean(VIPSOM_2P.CR(:,3:8),2);
VIPSOM_2P.CR_stages(:,3) = nanmean(VIPSOM_2P.CR(:,9:16),2);
VIPSOM_2P.CR_stages(:,4) = nanmean(VIPSOM_2P.CR(:,17:22),2);
temp_var = VIPSOM_2P.RwdMov_Duration;
temp_var(VIPSOM_2P.TrialNum<trialnum_thre) = nan;
VIPSOM_2P.RwdMov_Duration_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.RwdMov_Duration_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.RwdMov_Duration_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.RwdMov_Duration_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.RT;
VIPSOM_2P.RT_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.RT_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.RT_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.RT_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.RT_var;
VIPSOM_2P.RT_var_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.RT_var_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.RT_var_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.RT_var_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.RwdMVM_Rwd;
temp_var(VIPSOM_2P.TrialNum<trialnum_thre) = nan;
VIPSOM_2P.RwdMVM_Rwd_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.RwdMVM_Rwd_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.RwdMVM_Rwd_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.RwdMVM_Rwd_stages(:,4) = nanmean(temp_var(:,17:22),2);
for ii = 1:size(VIPSOM_2P.LeverCorr_Reward_Matrix,3)
    VIPSOM_2P.LeverCorr_Reward_within(ii,:) = diag(VIPSOM_2P.LeverCorr_Reward_Matrix(:,:,ii),0);
    VIPSOM_2P.LeverCorr_Reward_across(ii,:) = diag(VIPSOM_2P.LeverCorr_Reward_Matrix(:,:,ii),1);
    VIPSOM_2P.LeverCorr_First_within(ii,:) = diag(VIPSOM_2P.LeverCorr_First_Matrix(:,:,ii),0);
    VIPSOM_2P.LeverCorr_First_across(ii,:) = diag(VIPSOM_2P.LeverCorr_First_Matrix(:,:,ii),1);
    VIPSOM_2P.LeverCorr_All_within(ii,:) = diag(VIPSOM_2P.LeverCorr_All_Matrix(:,:,ii),0);
    VIPSOM_2P.LeverCorr_All_across(ii,:) = diag(VIPSOM_2P.LeverCorr_All_Matrix(:,:,ii),1);
end
temp_var = VIPSOM_2P.LeverCorr_Reward_within;
temp_var(VIPSOM_2P.TrialNum<trialnum_thre) = nan;
VIPSOM_2P.LeverCorr_within_Reward_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.LeverCorr_within_Reward_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.LeverCorr_within_Reward_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.LeverCorr_within_Reward_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.LeverCorr_Reward_across;
VIPSOM_2P.LeverCorr_across_Reward_stages(:,1) = nanmean(temp_var(:,1),2);
VIPSOM_2P.LeverCorr_across_Reward_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIPSOM_2P.LeverCorr_across_Reward_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIPSOM_2P.LeverCorr_across_Reward_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = VIPSOM_2P.LeverCorr_First_within;
VIPSOM_2P.LeverCorr_within_First_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.LeverCorr_within_First_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.LeverCorr_within_First_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.LeverCorr_within_First_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.LeverCorr_First_across;
VIPSOM_2P.LeverCorr_across_First_stages(:,1) = nanmean(temp_var(:,1),2);
VIPSOM_2P.LeverCorr_across_First_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIPSOM_2P.LeverCorr_across_First_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIPSOM_2P.LeverCorr_across_First_stages(:,4) = nanmean(temp_var(:,17:21),2);
temp_var = VIPSOM_2P.LeverCorr_All_within;
VIPSOM_2P.LeverCorr_within_All_stages(:,1) = nanmean(temp_var(:,1:2),2);
VIPSOM_2P.LeverCorr_within_All_stages(:,2) = nanmean(temp_var(:,3:8),2);
VIPSOM_2P.LeverCorr_within_All_stages(:,3) = nanmean(temp_var(:,9:16),2);
VIPSOM_2P.LeverCorr_within_All_stages(:,4) = nanmean(temp_var(:,17:22),2);
temp_var = VIPSOM_2P.LeverCorr_All_across;
VIPSOM_2P.LeverCorr_across_All_stages(:,1) = nanmean(temp_var(:,1),2);
VIPSOM_2P.LeverCorr_across_All_stages(:,2) = nanmean(temp_var(:,3:7),2);
VIPSOM_2P.LeverCorr_across_All_stages(:,3) = nanmean(temp_var(:,9:15),2);
VIPSOM_2P.LeverCorr_across_All_stages(:,4) = nanmean(temp_var(:,17:21),2);

save(['Z:\People\Chi\TwoP_IN\BhaviorAnalysis' filesep 'BehaviorAnalysis'],'SOM_2P','VIP_2P','-v7.3');

