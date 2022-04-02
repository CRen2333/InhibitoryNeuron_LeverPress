%% Collect_behavior_2P
clear;
close all;
clc;

IN = 'VIP_SOM';
Animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R',...
    'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
Stage = '';

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage '\MovAnalysis\' Animal '_CuedMov_All.mat']);
    
    for curr_day = 3:4 % collect day 3 and day 4 data
        ActiveTrials_binary{curr_animal,curr_day} = [];
        RwdTrials_binary{curr_animal,curr_day} = [];
        RspdTrials_binary{curr_animal,curr_day} = [];
        RwdMvMTrials_binary{curr_animal,curr_day} = [];
        CuedRwdMvMTrials_binary{curr_animal,curr_day} = [];        
        RT{curr_animal,curr_day} = [];
        RwdMov_Duration{curr_animal,curr_day} = [];
        C_CRM{curr_animal,curr_day} = [];
        RwdMVM_Rwd{curr_animal,curr_day} = [];
        LeverCorrRwd_Matrix{curr_animal,curr_day} = [];
        RwdMVM_order{curr_animal,curr_day} = [];
        CuedRwdMVM_order{curr_animal,curr_day} = [];
        Attemp_order{curr_animal,curr_day} = [];
        total_trialnum = [];
        if isempty(CuedMov_SingleAnimal{curr_day})
            continue
        end
        % save trial index
        for curr_session = 1:length(CuedMov_SingleAnimal{curr_day}.Bpod)
            total_trialnum(curr_session) = length(CuedMov_SingleAnimal{curr_day}.Bpod{curr_session});
            temp_info = CuedMov_SingleAnimal{curr_day}.Cued_MovOnset_Info{curr_session};
            if curr_session > 1
                temp_info(:,1) = temp_info(:,1)+sum(total_trialnum(1:curr_session-1));
            end
            first_index = temp_info(:,4)==1;
            rewarded_mvm_index = ~isnan(temp_info(:,7));
            cued_rewarded_mvm_index = temp_info(:,5)==1;
        
            rspd_trial = temp_info(first_index,1);
            rwdmvm_trial = temp_info(rewarded_mvm_index,1);
            cued_rwdmvm_trial = temp_info(cued_rewarded_mvm_index,1);
            
            if curr_session > 1
                rwd_trial_last = rwd_trial(end);
                rwd_trial = find(~isnan(CuedMov_SingleAnimal{curr_day}.Bpod{curr_session}(:,5)));
                rwd_trial = rwd_trial + sum(total_trialnum(1:curr_session-1));
                cued_rwdmvm_in_rwd = find(ismember(rwd_trial,cued_rwdmvm_trial))+rwd_trial_last;
            else
                rwd_trial = find(~isnan(CuedMov_SingleAnimal{curr_day}.Bpod{curr_session}(:,5)));
                cued_rwdmvm_in_rwd = find(ismember(rwd_trial,cued_rwdmvm_trial));
            end
            
            if curr_session > 1
                RspdTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RspdTrials_binary{curr_animal,curr_day}{curr_session,1}(rspd_trial-sum(total_trialnum(1:curr_session-1))) = 1;
                RwdTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RwdTrials_binary{curr_animal,curr_day}{curr_session,1}(rwd_trial-sum(total_trialnum(1:curr_session-1))) = 1;
                RwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1}(rwdmvm_trial-sum(total_trialnum(1:curr_session-1))) = 1;
                CuedRwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                CuedRwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1}(cued_rwdmvm_trial-sum(total_trialnum(1:curr_session-1))) = 1;
                ActiveTrials_binary{curr_animal,curr_day}{curr_session,1} = logical(RspdTrials_binary{curr_animal,curr_day}{curr_session,1}+RwdTrials_binary{curr_animal,curr_day}{curr_session,1});
            else
                RspdTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RspdTrials_binary{curr_animal,curr_day}{curr_session,1}(rspd_trial) = 1;
                RwdTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RwdTrials_binary{curr_animal,curr_day}{curr_session,1}(rwd_trial) = 1;
                RwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                RwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1}(rwdmvm_trial) = 1;
                CuedRwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1} = zeros(total_trialnum(curr_session),1);
                CuedRwdMvMTrials_binary{curr_animal,curr_day}{curr_session,1}(cued_rwdmvm_trial) = 1;
                ActiveTrials_binary{curr_animal,curr_day}{curr_session,1} = logical(RspdTrials_binary{curr_animal,curr_day}{curr_session,1}+RwdTrials_binary{curr_animal,curr_day}{curr_session,1});
            end
           

            % Trace for moving average
            
            % reaction time
            RT{curr_animal,curr_day}{curr_session,1}(:,1) = rspd_trial;
            RT{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(first_index,2)-temp_info(first_index,6);
            RT{curr_animal,curr_day}{curr_session,1}(:,3) = temp_info(first_index,5);
            RT{curr_animal,curr_day}{curr_session,1}(:,4) = ~isnan(temp_info(first_index,7));
            % rewarded movement related
            RwdMov_Duration{curr_animal,curr_day}{curr_session,1}(:,1) = cued_rwdmvm_trial;
            RwdMov_Duration{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(cued_rewarded_mvm_index,3);
            RwdMov_Duration{curr_animal,curr_day}{curr_session,1}(:,3) = cued_rwdmvm_in_rwd;
            C_CRM{curr_animal,curr_day}{curr_session,1}(:,1) = cued_rwdmvm_trial;
            C_CRM{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(cued_rewarded_mvm_index,2)-temp_info(cued_rewarded_mvm_index,6);
            C_CRM{curr_animal,curr_day}{curr_session,1}(:,3) = cued_rwdmvm_in_rwd;
            nan_index = (temp_info(cued_rewarded_mvm_index,2)-temp_info(cued_rewarded_mvm_index,6))>10;
            C_CRM{curr_animal,curr_day}{curr_session,1}(nan_index,2) = nan;
            RwdMVM_Rwd{curr_animal,curr_day}{curr_session,1}(:,1) = cued_rwdmvm_trial;
            RwdMVM_Rwd{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(logical(cued_rewarded_mvm_index),7)-temp_info(logical(cued_rewarded_mvm_index),2);
            RwdMVM_Rwd{curr_animal,curr_day}{curr_session,1}(:,3) = cued_rwdmvm_in_rwd;
            nan_index = temp_info(logical(cued_rewarded_mvm_index),7)-temp_info(logical(cued_rewarded_mvm_index),2)>10;
            RwdMVM_Rwd{curr_animal,curr_day}{curr_session,1}(nan_index,2) = nan;
            RwdMVM_order{curr_animal,curr_day}{curr_session,1}(:,1) = rwdmvm_trial; % trial index
            RwdMVM_order{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(rewarded_mvm_index,4); % mvm order
            CuedRwdMVM_order{curr_animal,curr_day}{curr_session,1}(:,1) = cued_rwdmvm_trial; % trial index
            CuedRwdMVM_order{curr_animal,curr_day}{curr_session,1}(:,2) = temp_info(cued_rewarded_mvm_index,4); % mvm order
            CuedRwdMVM_order{curr_animal,curr_day}{curr_session,1}(:,3) = cued_rwdmvm_in_rwd;
            trial_index = unique(temp_info(:,1));
            Attemp_order{curr_animal,curr_day}{curr_session,1} = [];
            for trial_ii = 1:length(trial_index)
                curr_index = find(temp_info(:,1)==trial_index(trial_ii),1,'last');
                Attemp_order{curr_animal,curr_day}{curr_session,1} = [Attemp_order{curr_animal,curr_day}{curr_session,1};temp_info(curr_index,[1,4,5,7])];
            end
                           
        end
        RspdTrials_binary{curr_animal,curr_day} = cell2mat(RspdTrials_binary{curr_animal,curr_day});
        RwdTrials_binary{curr_animal,curr_day} = cell2mat(RwdTrials_binary{curr_animal,curr_day});
        RwdMvMTrials_binary{curr_animal,curr_day} = cell2mat(RwdMvMTrials_binary{curr_animal,curr_day});
        CuedRwdMvMTrials_binary{curr_animal,curr_day} = cell2mat(CuedRwdMvMTrials_binary{curr_animal,curr_day});
        ActiveTrials_binary{curr_animal,curr_day} = cell2mat(ActiveTrials_binary{curr_animal,curr_day});
        RT{curr_animal,curr_day} = cell2mat(RT{curr_animal,curr_day});
        RwdMov_Duration{curr_animal,curr_day} = cell2mat(RwdMov_Duration{curr_animal,curr_day});
        C_CRM{curr_animal,curr_day} = cell2mat(C_CRM{curr_animal,curr_day});
        RwdMVM_Rwd{curr_animal,curr_day} = cell2mat(RwdMVM_Rwd{curr_animal,curr_day});
        RwdMVM_order{curr_animal,curr_day} = cell2mat(RwdMVM_order{curr_animal,curr_day});
        CuedRwdMVM_order{curr_animal,curr_day} = cell2mat(CuedRwdMVM_order{curr_animal,curr_day});
        Attemp_order{curr_animal,curr_day} = cell2mat(Attemp_order{curr_animal,curr_day});
        
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
        LeverCorrRwd_Matrix{curr_animal,curr_day} = corrcoef(temp_trace);
                        
    end
        
    clear temp_trace Trial_Trial_Corr_Reward Trial_Trial_Corr_First Trial_Trial_Corr_All CuedMov_SingleAnimal
end

save('BehaviorCollect_Nai_Learning_Day3.mat','-v7.3');

%% hM4Di
Animals = Animals(1:22);
hM4Di_animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R'};
mCherry_animals = {'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
hM4Di_index = ismember(Animals,hM4Di_animals);
mCherry_index = ismember(Animals,mCherry_animals);

% CR, first 10 trials
temp_var = [];
for curr_animal = 1:length(Animals)
    temp_var(curr_animal,1) = sum(RwdTrials_binary{curr_animal,3}(1:10))/10;
end

% RwdMVM_Rwd, first 10 trials
temp_var = [];
for curr_animal = 1:length(Animals)
    index = RwdMVM_Rwd{curr_animal,3}(:,1)<=10; % first 10 traisl
    temp_var(curr_animal,1) = nanmedian(RwdMVM_Rwd{curr_animal,3}(index,2));
end

temp_var_1 = temp_var(mCherry_index,:);
temp_var_2 = temp_var(hM4Di_index,:);
[CX] = cbrewer('seq','RdPu',9);
color_value = nanmean(CX(6:7,:));    

figure; hold on; set(gcf,'position',[200 200 200 200],'color','w');
plot(1+((rand(1,11)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,11)-0.5)*0.25),temp_var_2(:),'color',color_value,'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor',color_value,'FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color',[0.5 0.5 0.5],'linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color',color_value,'linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'mCherry','hM4Di'});
% ylim([0 1]);

kk = 1;
Drug = repmat(hM4Di_index',1,1);
Drug = Drug(:);
Sessions = repmat([1:kk],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,kk);
Animals_test = Animals_test(:);
y = temp_var(:,1:kk);
y = y(:);
tbl = table(Drug,Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
line([0.9,2.1],[1 1],'color','k')
text(1.5,1,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

% ylim([0 1]); ylabel('Correct rate');
% ylim([0 3]); ylabel('RwdMVM.Onset to Rwd (s)');
axis square



