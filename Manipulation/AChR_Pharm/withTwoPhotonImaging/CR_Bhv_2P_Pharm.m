%% Collect_behavior_2P
clearvars -except SOM VIP;
close all;
clc;

IN = 'SOM'; 
SessionType = 'Exp_Ago';
switch IN
    case 'VIP'
        switch SessionType
            case 'Nai_Ant'
                Animals = {'CR_3702608-O','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3887044-L','CR_3786160-LR','CR_3887043-O','CR_3887043-L','CR_3887041-O'};
                Drug = 'Ant';
            case 'Exp_Ago'
                Animals = {'WL_3526642-L','WL_3526642-R','CR_3619106-R','CR_3619106-LR','CR_3633192-L','CR_3633192-R','CR_3702608-O','CR_3702608-LR','CR_3658844-R','CR_3672031-O','CR_3633193-L'};           
                Drug = 'Ago';
        end        
    case 'SOM'
        switch SessionType
            case 'Nai_Ant'
                Animals = {'CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O','CR_3886982-O','CR_3886982-L','CR_3886982-LR','CR_3886982-R','CR_3887040-O','CR_3887040-L'};
                Drug = 'Ant';
            case 'Exp_Ago'
                Animals = {'WL_3547273-LR','WL_3547272-O','WL_3547272-L','WL_3526578-O','WL_3526578-R','CR_3672035-R','CR_3702224-L','CR_3702224-O','CR_3702224-R','CR_3702226-O'};
                Drug = 'Ago';
        end
end
       
Fields = {'Field_1','Field_2'};

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    switch SessionType
        case 'Exp_Ago'
            if ismember(Animal,{'CR_3633193-L','CR_3658844-L','CR_3658844-R','CR_3672031-O','CR_3672035-R','CR_3702224-L','CR_3702224-O','CR_3702224-R','CR_3702226-O'})
                SessionType_2 = 'Exp_Ago_Ant';
            else
                SessionType_2 = SessionType;
            end
        otherwise
            SessionType_2 = SessionType;
    end
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep SessionType_2 '\MovAnalysis\' Animal '_CuedMov_All.mat']);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep SessionType_2 '\df_f\' Animal '_ImagingInfo.mat'],'Imaging_Fields');

    Days = length(Imaging_Fields{1}.Date);
    RwdMov_Duration(curr_animal,[11:2]) = nan;
    CR(curr_animal,[1:2]) = nan;
    RT(curr_animal,[1:2]) = nan;
    C_CRM(curr_animal,[1:2]) = nan;
    RwdMVM_Rwd(curr_animal,[1:2]) = nan;
    TrialNum(curr_animal,[1:2]) = nan;
    Bhv_dates = dir(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep SessionType_2 filesep 'LeverTrace']);
    Bhv_dates = {Bhv_dates.name};
    Bhv_dates = Bhv_dates(3:end);
    imaging_index = ismember(Bhv_dates,Imaging_Fields{1}.Date);
    CuedMov_SingleAnimal = CuedMov_SingleAnimal(imaging_index);
    Trial_Trial_Corr_Reward = Trial_Trial_Corr_Reward(imaging_index,imaging_index);
    for curr_day = 1:2
        if isempty(CuedMov_SingleAnimal{curr_day})
            continue
        end
        RwdMov_Duration(curr_animal,curr_day) = CuedMov_SingleAnimal{curr_day}.Median_Movduration_Reward;
        if strcmp(SessionType,'Exp_Ago') && ismember(Animal,{'CR_3619106-R','CR_3619106-LR','CR_3633192-L','CR_3633192-R','CR_3633193-L'})
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
        C_CRM(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,2)-temp_info(rewared_mvm_index,6));
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
        TrialNum(curr_animal,curr_day) = size(temp_trace,2);
        RwdMVM_Rwd(curr_animal,curr_day) = nanmedian(temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2));
        clear temp_trace
    end
    Corr_within(curr_animal,:) = diag(Trial_Trial_Corr_Reward([1,2],[1,2]));
    clear Trial_Trial_Corr_Reward Trial_Trial_Corr_First Trial_Trial_Corr_All CuedMov_SingleAnimal
    if strcmp(SessionType,'Nai_Ant')
        if strcmp(Imaging_Fields{1}.SessionType{1},'Nai_Ant')
            RwdMov_Duration(curr_animal,:) = fliplr(RwdMov_Duration(curr_animal,:));
            CR(curr_animal,:) = fliplr(CR(curr_animal,:));
            RT(curr_animal,:) = fliplr(RT(curr_animal,:));
            C_CRM(curr_animal,:) = fliplr(C_CRM(curr_animal,:));
            RwdMVM_Rwd(curr_animal,:) = fliplr(RwdMVM_Rwd(curr_animal,:));
            TrialNum(curr_animal,:) = fliplr(TrialNum(curr_animal,:));
            Corr_within(curr_animal,:) = fliplr(Corr_within(curr_animal,:));
        end
    end
    clear Imajing_Fields
end

switch IN
    case 'SOM'
        SOM.RwdMov_Duration = RwdMov_Duration;
        SOM.CR = CR;
        SOM.RT = RT;
        SOM.RwdMVM_Rwd = RwdMVM_Rwd;
        SOM.Corr_within = Corr_within;
        SOM.C_CRM = C_CRM;
        SOM.TrialNum = TrialNum;
    case 'VIP'
        VIP.RwdMov_Duration = RwdMov_Duration;
        VIP.CR = CR;
        VIP.RT = RT;
        VIP.RwdMVM_Rwd = RwdMVM_Rwd;
        VIP.Corr_within = Corr_within;
        VIP.C_CRM = C_CRM;
        VIP.TrialNum = TrialNum;
        
end
save(['Z:\People\Chi\TwoP_IN\BhaviorAnalysis' filesep 'BehaviorAnalysis_NaiAnt'],'SOM','VIP','-v7.3');
% save(['Z:\People\Chi\TwoP_IN\BhaviorAnalysis' filesep 'BehaviorAnalysis_ExpAgo'],'SOM','VIP','-v7.3');

%% Plot, combine VIP and SOM
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.Corr_within;SOM.Corr_within];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');

% p = friedman([temp_var_1(:),temp_var_2(:)],6)
% ranksum(temp_var_1(:,4),temp_var_2(:,4))
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.CR = pValue(2);
text(1.5,1.2,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');
axis square;

% ylabel('Correct rate'); ylim([0 1]);
% ylabel('Reaction time (s)'); 
% ylabel('PC1 Fraction'); 
% ylabel('Mvm. duration (sec)'); ylim([1 3.5]);
% ylabel({'C_CuedRwdMvmOnset (s)'});
% ylabel({'std. CuedRwdMvmOnset', 'to Rwd (s)'});
% ylabel('Mvm. speed (abs.)');
% ylabel('std. Reaction time (s)'); 
% ylabel({'std. Cue to', 'CuedRwdMvmOnset (s)'});
% ylabel('Corr. coef (across)'); ylim([0 0.3])
% ylabel('Lever corr. within (rwd)'); ylim([0 0.35])
ylabel('CuedRwdMvmOnset_R (s)'); ylim([0 1.2])
% ylabel({'CuedRwdMvmOnset','to Reward (s)'});
% ylabel({'std. CuedRwdMvmOnset','to Reward (s)'});
% ylabel('std. # of movements');
% ylabel('# of movements (median)');
% ylabel('# of movements (mean)');

