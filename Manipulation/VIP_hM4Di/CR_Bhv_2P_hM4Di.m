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
    temp_matrix = nan(2,2);
    temp_matrix(1:size(Trial_Trial_Corr_All,1),1:size(Trial_Trial_Corr_All,1)) = Trial_Trial_Corr_All;
    LeverCorr_All_Matrix(:,:,curr_animal) = temp_matrix;
    temp_matrix = nan(2,2);
    temp_matrix(1:size(Trial_Trial_Corr_First,1),1:size(Trial_Trial_Corr_First,1)) = Trial_Trial_Corr_First;
    LeverCorr_First_Matrix(:,:,curr_animal) = temp_matrix;
    RwdMov_Duration(curr_animal,[1:2]) = nan;
    CR(curr_animal,[1:2]) = nan;
    RT(curr_animal,[1:2]) = nan;
    RT_var(curr_animal,[1:2]) = nan;
    RwdMVM_Rwd(curr_animal,[1:2]) = nan;
    RwdMVM_order(curr_animal,[1:2]) = nan;
    RwdMVM_order_mean(curr_animal,[1:2]) = nan;
    C_CRM(curr_animal,[1:2]) = nan;
    First_R(curr_animal,[1:2]) = nan;
    TrialNum(curr_animal,[1:2]) = nan;
    RMS_reward(curr_animal,[1:2]) = nan;
    RMS_all(curr_animal,[1:2]) = nan;
    RMS_first(curr_animal,[1:2]) = nan;
    RwdMVM_speed(curr_animal,[1:2]) = nan;
    FirstMVM_speed(curr_animal,[1:2]) = nan;
    AllMVM_speed(curr_animal,[1:2]) = nan;
    
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
        RT_var(curr_animal,curr_day) = nanstd(temp_info(first_index,2)-temp_info(first_index,6));
        RwdMVM_Rwd(curr_animal,curr_day) = nanmedian(temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2));
        C_CRM(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,2)-temp_info(rewared_mvm_index,6));
        First_R(curr_animal,curr_day) = nanmedian(temp_info(first_index,7)-temp_info(first_index,2));
        RwdMVM_order(curr_animal,curr_day) = nanmedian(temp_info(rewared_mvm_index,4));
        RwdMVM_order_mean(curr_animal,curr_day) = nanmean(temp_info(rewared_mvm_index,4));
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_Reward_downsample;
        temp_RwdMVM_Rwd = temp_info(logical(reward_index.*rewared_mvm_index),7)-temp_info(logical(reward_index.*rewared_mvm_index),2);
        TrialNum(curr_animal,curr_day) = size(temp_trace,2);
        temp_RwdMVM_Rwd = round(temp_RwdMVM_Rwd*100);
        for tt = 1:TrialNum(curr_animal,curr_day)
            temp_RwdMVM_speed(tt) = nanmean(abs(diff(temp_trace(1:min(200,temp_RwdMVM_Rwd(tt)),tt))));
        end
        RwdMVM_speed(curr_animal,curr_day) = nanmean(temp_RwdMVM_speed);
        clear temp_RwdMVM_Rwd temp_RwdMVM_speed        

        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_All_downsample;
        AllMVM_speed(curr_animal,curr_day) = nanmean(nanmean(abs(diff(temp_trace,[],1)),1))/10;
        temp_trace = CuedMov_SingleAnimal{curr_day}.LeverTrace_First_downsample;
        FirstMVM_speed(curr_animal,curr_day) = nanmean(nanmean(abs(diff(temp_trace,[],1)),1))/10;
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

for ii = 1:length(Animals)
    VIPSOM_2P.LeverCorr_Reward_within(ii,:) = diag(LeverCorr_Reward_Matrix(:,:,ii));
end

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
VIPSOM_2P.RwdMVM_order = RwdMVM_order;
VIPSOM_2P.RwdMVM_order = RwdMVM_order;
VIPSOM_2P.TrialNum = TrialNum;
VIPSOM_2P.RMS_reward = RMS_reward;
VIPSOM_2P.RMS_all = RMS_all;
VIPSOM_2P.RMS_first = RMS_first;
VIPSOM_2P.RwdMVM_speed = RwdMVM_speed;
VIPSOM_2P.AllMVM_speed = AllMVM_speed;
VIPSOM_2P.FirstMVM_speed = FirstMVM_speed;

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
%             temp_lever_trace_session = resample(Lever_force_resample{jj},1,10);
%             temp_lever_trace_session([1:10,end-9:end]) = nan;
            butterworth_stop = 5/500; % fraction of nyquist (cutoff = 10 Hz)
            [b a] = butter(4, butterworth_stop,'low');
            lever_force_smooth = filtfilt(b,a,Lever_force_resample{jj});
            lever_velocity_resample = [0;diff(lever_force_smooth)];
            lever_velocity_resample_smooth = smooth(lever_velocity_resample,1);
            time_info = round(time_info*1000);            
            for tt = 1:length(time_info)
%                 speed_mvm_rwd{jj,1}(tt,1) = nanmean(abs(diff((temp_lever_trace_session(time_info(tt,1):time_info(tt,2))))));
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
        first_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,4) == 1;
        rwd_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5) == 1;
        temp_matrix = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(first_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(first_index,2);
        First_Rwd(curr_animal,ii) = nanmedian(temp_matrix);
        clear temp_matrix
        hM4Di_rwd_in_first(curr_animal,ii) = sum(first_index.*rwd_index)/sum(first_index);
        % Traces
        hM4Di_Trace_all{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_All,1,10);
        hM4Di_Trace_all{curr_animal,ii} = hM4Di_Trace_all{curr_animal,ii}(10:end-9,:);
        hM4Di_Trace_first{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_First,1,10);
        hM4Di_Trace_first{curr_animal,ii} = hM4Di_Trace_first{curr_animal,ii}(10:end-9,:);
        hM4Di_Trace_reward{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_Reward,1,10);
        hM4Di_Trace_reward{curr_animal,ii} = hM4Di_Trace_reward{curr_animal,ii}(10:end-9,:);
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
%             temp_lever_trace_session = resample(Lever_force_resample{jj},1,10);
%             temp_lever_trace_session([1:10,end-9:end]) = nan;
            butterworth_stop = 5/500; % fraction of nyquist (cutoff = 10 Hz)
            [b a] = butter(4, butterworth_stop,'low');
            lever_force_smooth = filtfilt(b,a,Lever_force_resample{jj});
            lever_velocity_resample = [0;diff(lever_force_smooth)];
            lever_velocity_resample_smooth = smooth(lever_velocity_resample,1);
            time_info = round(time_info*1000);            
            for tt = 1:length(time_info)
%                 speed_mvm_rwd{jj,1}(tt,1) = nanmean(abs(diff((temp_lever_trace_session(time_info(tt,1):time_info(tt,2))))));
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
        first_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,4) == 1;
        rwd_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5) == 1;
        mCherry_rwd_in_first(curr_animal,ii) = sum(first_index.*rwd_index)/sum(first_index);
        % Traces
        mCherry_Trace_all{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_All,1,10);
        mCherry_Trace_all{curr_animal,ii} = mCherry_Trace_all{curr_animal,ii}(10:end-9,:);
        mCherry_Trace_first{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_First,1,10);
        mCherry_Trace_first{curr_animal,ii} = mCherry_Trace_first{curr_animal,ii}(10:end-9,:);
        mCherry_Trace_reward{curr_animal,ii} = resample(CuedMov_SingleAnimal{ii}.LeverTrace_Reward,1,10);
        mCherry_Trace_reward{curr_animal,ii} = mCherry_Trace_reward{curr_animal,ii}(10:end-9,:);
    end
    clear CuedMov_SingleAnimal
end

%% hM4Di
temp_var = VIPSOM_2P.CR;
temp_var_1 = temp_var(mCherry_index,1:2);
temp_var_2 = temp_var(hM4Di_index,1:2);
[CX] = cbrewer('seq','RdPu',9);
color_value = nanmean(CX(6:7,:));    

figure; hold on; set(gcf,'position',[200 200 200 200],'color','w');
plot(1+((rand(1,22)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,22)-0.5)*0.25),temp_var_2(:),'color',color_value,'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor',color_value,'FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color',[0.5 0.5 0.5],'linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color',color_value,'linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'mCherry','hM4Di'});
ylim([0 9]); yticks([0:3:9]);

% ylim([0 1]); ylabel('Correct rate');
% ylim([0 3.3]); ylabel('Reaction time (s)');
ylim([0 6]); ylabel('Cue to RwdMVM.Onset (s)');
% ylim([0 2.5]); ylabel('RwdMVM.Onset to Rwd (s)');
% ylim([0 8]); ylabel('RwdMVM Duration (s)');
% ylim([0 0.55]); ylabel('Lever Corr. within (rwd)');
ylim([0.25 0.6]); ylabel('PC1 fraction (trial as var)');
ylim([0 0.085]); ylabel('Rwd MVM speed (rwd)');
ylabel('1st Mvm to Rwd (s)')
axis square

Drug = repmat(hM4Di_index',1,2);
Drug = Drug(:);
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:,1:2);
y = y(:);
tbl = table(Drug,Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
line([0.5,1.5],[0.4 0.4],'color','k')
text(1,0.42,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);

% movement speed       
temp_var_1 = mCherry_speed_mvm_rwd*1000*7.2;
temp_var_2 = hM4Di_speed_mvm_rwd*1000*7.2;
figure; hold on; set(gcf,'position',[200 200 200 200],'color','w');
plot(1+((rand(1,22)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,22)-0.5)*0.25),temp_var_2(:),'color',color_value,'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor',color_value,'FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color',[0.5 0.5 0.5],'linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color',color_value,'linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'mCherry','hM4Di'});
ylim([0 15]); 
axis square;

ylabel('Fraction of rspd. trials')
ylabel('Fraction of active epoch')
ylabel('Fraction of fst. rwd./rwd.')
ylabel('Speed (mm/s)');
temp_var = [temp_var_1;temp_var_2];
Drug = repmat(~hM4Di_index',1,2);
Drug = Drug(:);
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:,1:2);
y = y(:);
tbl = table(Drug,Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,15,['hM4Di: ' num2str(pValue(2))],'fontsize',6);

figure; set(gcf,'color','w','position',[200,200,800,800])
hold on;
for ii = 1:size(mCherry_Trace_all,1)
    for jj = 1:size(mCherry_Trace_all,2)
        temp_trace = mCherry_Trace_reward{ii,jj}(1:240,:);
        subplot(8,8,jj+(ii-1)*8); hold on;
        plot(temp_trace,'color',[0.7,0.7,0.7],'linewidth',0.2);
        plot(nanmean(temp_trace,2),'k','linewidth',1);
        ymin = min(nanmean(temp_trace,2))-0.5;
        ymax = max(nanmean(temp_trace,2))+0.5;
        xlim([1 240]); ylim([ymin ymax]),line([40 40],ylim,'color','k','linestyle',':');
        axis off
    end    
end

figure; set(gcf,'color','w','position',[200,200,800,800])
hold on;
for ii = 1:size(hM4Di_Trace_all,1)
    for jj = 1:size(hM4Di_Trace_all,2)
        temp_trace = hM4Di_Trace_all{ii,jj}(1:240,:);
        subplot(8,8,jj+(ii-1)*8); hold on;
        plot(temp_trace,'color',[0.7,0.7,0.7],'linewidth',0.2);
        plot(nanmean(temp_trace,2),'r','linewidth',1);
        ymin = min(nanmean(temp_trace,2))-0.5;
        ymax = max(nanmean(temp_trace,2))+0.5;
        xlim([1 240]); ylim([ymin ymax]),line([40 40],ylim,'color','k','linestyle',':');
        axis off
    end    
end

% Correlation matrix
% Reward
temp_var = VIPSOM_2P.LeverCorr_Reward_Matrix;
for ii = 1:size(temp_var,3)
    temp_trial_num = VIPSOM_2P.TrialNum(ii,:)<trialnum_thre;
    temp_trial_num = repmat(temp_trial_num,length(temp_trial_num),1)+repmat(temp_trial_num',1,length(temp_trial_num));
    temp_trial_num = logical(temp_trial_num);
    temp_temp_var = temp_var(:,:,ii);
    temp_temp_var(temp_trial_num) = nan;
    temp_var(:,:,ii) = temp_temp_var;
end

temp_var_1 = temp_var([1:8],[1:8],hM4Di_index);
temp_var_2 = temp_var([1:8],[1:8],mCherry_index);

hM4Di_plot = temp_var_1([1:2],:,:);
mCherry_plot = temp_var_2([1:2],:,:);

figure; set(gcf,'color','w','position',[200 200 600 200]);
hold on;
subplot(1,2,1);imagesc(nanmean(temp_var_1,3),[0,0.3]);colormap hot; xlabel('session');ylabel('session');
colorbar; title('hM4Di');
axis square;
subplot(1,2,2);imagesc(nanmean(temp_var_2,3),[0,0.3]);colormap hot; xlabel('session');ylabel('session');
axis square; title('mCherry');
colorbar

figure; set(gcf,'color','w','position',[200 200 300 200]);
hold on; c_color = cbrewer('div','RdBu',64);c_color = flipud(c_color);
imagesc(nanmean(temp_var_1,3)-nanmean(temp_var_2,3),[-0.15 0.15]);colormap(c_color); xlabel('session');ylabel('session');
colorbar;
axis square; xlim([0.5,8.5]);ylim([0.5 8.5]);set(gca,'YDir','reverse');box on;
title('hM4Di - mCherry')

test_session = 2;
figure; hold on; set(gcf,'position',[200 200 200 200],'color','w');
temp_temp_var_1 = squeeze(hM4Di_plot(test_session,:,:))';
temp_temp_var_2 = squeeze(mCherry_plot(test_session,:,:))';
plot(temp_temp_var_1','color',[1,0.7,0.7]);
plot(temp_temp_var_2','color',[0.7,0.7,0.7]);
temp_var_11 = nanmean(temp_temp_var_1,1);
temp_var_12 = nanstd(temp_temp_var_1,[],1)./sqrt(sum(~isnan(temp_temp_var_1)));
temp_var_21 = nanmean(temp_temp_var_2,1);
temp_var_22 = nanstd(temp_temp_var_2,[],1)./sqrt(sum(~isnan(temp_temp_var_2)));
plot(temp_var_11,'color','r','linewidth',1);
plot(temp_var_21,'color','k','linewidth',1);
for ii = 1:7
    line([ii,ii],[temp_var_11(ii)-temp_var_12(ii),temp_var_11(ii)+temp_var_12(ii)],'color','r','linewidth',1);
    line([ii,ii],[temp_var_21(ii)-temp_var_22(ii),temp_var_21(ii)+temp_var_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 8.5]); xlabel('session'); axis square; ylabel('Lever Corr. with session 2');

Drug = repmat(hM4Di_index',1,8);
Drug = Drug(:);
Sessions = repmat([1:8],size(Animals,2),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(Animals,2)]',1,8);
Animals_test = Animals_test(:);
y = temp_var(1:8,1:8,:);
y = y(test_session,:,:);
y = squeeze(y);
y = y';
y = y(:);
tbl = table(Drug,Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
line([1,8],[0.4,0.4],'color','k')
text(5,0.42,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);


