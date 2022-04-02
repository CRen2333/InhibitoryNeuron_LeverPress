%% combine bhv from 2P and WF
clear
load('Z:\People\Chi\WFLP_IN\All_BehaviorCollect_withPupilAnimal_refined_ForBhv.mat','VIP','SOM','-mat');
load('Z:\People\Chi\TwoP_IN\BhaviorAnalysis\BehaviorAnalysis.mat');

% combine
session_control = [1:22];
session_control_2 = [1:21];
SOM_Combine.CR = [SOM.CR(:,session_control);SOM_2P.CR];
SOM_Combine.RT = [SOM.Reaction_Time(:,session_control);SOM_2P.RT];
SOM_Combine.Cue_to_CuedRewardedMov = [SOM.Cue_to_CuedRewardedMov(:,session_control);SOM_2P.C_CRM];
SOM_Combine.Cue_to_Reward = [SOM.Cue_to_Reward(:,session_control);SOM_2P.C_Rwd];
SOM_Combine.RwdMov_Duration = [SOM.CRM_Dur(:,session_control);SOM_2P.RwdMov_Duration];
SOM_Combine.LeverCorr_Reward_within = [SOM.Corr_Reward_within(:,session_control);SOM_2P.LeverCorr_Reward_within];
SOM_Combine.LeverCorr_Reward_across = [SOM.Corr_Reward_across(:,session_control_2);SOM_2P.LeverCorr_Reward_across];
SOM_Combine.RwdMVM_Rwd = [SOM.CuedRewardedMov_to_Reward(:,session_control);SOM_2P.RwdMVM_Rwd];
SOM_Combine.MovSpeed = [SOM.MovSpeed(:,session_control);SOM_2P.RwdMVM_speed];
SOM_Combine.TrialNum = [SOM.TrialNum(:,session_control);SOM_2P.TrialNum];

VIP_Combine.CR = [VIP.CR(:,session_control);VIP_2P.CR];
VIP_Combine.RT = [VIP.Reaction_Time(:,session_control);VIP_2P.RT];
VIP_Combine.Cue_to_CuedRewardedMov = [VIP.Cue_to_CuedRewardedMov(:,session_control);VIP_2P.C_CRM];
VIP_Combine.Cue_to_Reward = [VIP.Cue_to_Reward(:,session_control);VIP_2P.C_Rwd];
VIP_Combine.RwdMov_Duration = [VIP.CRM_Dur(:,session_control);VIP_2P.RwdMov_Duration];
VIP_Combine.LeverCorr_Reward_within = [VIP.Corr_Reward_within(:,session_control);VIP_2P.LeverCorr_Reward_within];
VIP_Combine.LeverCorr_Reward_across = [VIP.Corr_Reward_across(:,session_control_2);VIP_2P.LeverCorr_Reward_across];
VIP_Combine.RwdMVM_Rwd = [VIP.CuedRewardedMov_to_Reward(:,session_control);VIP_2P.RwdMVM_Rwd];
VIP_Combine.MovSpeed = [VIP.MovSpeed(:,session_control);VIP_2P.RwdMVM_speed];
VIP_Combine.TrialNum = [VIP.TrialNum(:,session_control);VIP_2P.TrialNum];

Combine.CR = [VIP_Combine.CR;SOM_Combine.CR];
Combine.RT = [VIP_Combine.RT;SOM_Combine.RT];
Combine.Cue_to_CuedRewardedMov = [VIP_Combine.Cue_to_CuedRewardedMov;SOM_Combine.Cue_to_CuedRewardedMov];
Combine.Cue_to_Reward = [VIP_Combine.Cue_to_Reward;SOM_Combine.Cue_to_Reward];
Combine.RwdMov_Duration = [VIP_Combine.RwdMov_Duration;SOM_Combine.RwdMov_Duration];
Combine.LeverCorr_Reward_within = [VIP_Combine.LeverCorr_Reward_within;SOM_Combine.LeverCorr_Reward_within];
Combine.LeverCorr_Reward_across = [VIP_Combine.LeverCorr_Reward_across;SOM_Combine.LeverCorr_Reward_across];
Combine.RwdMVM_Rwd = [VIP_Combine.RwdMVM_Rwd;SOM_Combine.RwdMVM_Rwd];
Combine.MovSpeed = [VIP_Combine.MovSpeed;SOM_Combine.MovSpeed];
Combine.TrialNum = [VIP_Combine.TrialNum;SOM_Combine.TrialNum];

% Matrix
SOM_Combine.LeverCorr_Reward_Matrix = [];
for ii = 1:size(SOM.Corr_Reward_Matrix,3)
    temp_matrix = SOM.Corr_Reward_Matrix(1:22,1:22,ii);
    temp_trialnum = SOM.TrialNum(ii,1:22);
    for m = 1:22
        for n = 1:22
            if temp_trialnum(m)<3 || temp_trialnum(n)<3
                temp_index(m,n) = true;
            else
                temp_index(m,n) = false;
            end
        end
    end
    temp_matrix(temp_index) = nan;
    SOM_Combine.LeverCorr_Reward_Matrix(:,:,ii) = temp_matrix;
    clear temp_matrix temp_trialnum temp_index
end
for jj = 1:size(SOM_2P.LeverCorr_Reward_Matrix,3)
    temp_matrix = SOM_2P.LeverCorr_Reward_Matrix(1:22,1:22,jj);
    temp_trialnum = SOM_2P.TrialNum(jj,1:22);
    for m = 1:22
        for n = 1:22
            if temp_trialnum(m)<3 || temp_trialnum(n)<3
                temp_index(m,n) = true;
            else
                temp_index(m,n) = false;
            end
        end
    end
    temp_matrix(temp_index) = nan;
    SOM_Combine.LeverCorr_Reward_Matrix(:,:,ii+jj) = temp_matrix;
    clear temp_matrix temp_trialnum temp_index
end

VIP_Combine.LeverCorr_Reward_Matrix = [];
for ii = 1:size(VIP.Corr_Reward_Matrix,3)
    temp_matrix = VIP.Corr_Reward_Matrix(1:22,1:22,ii);
    temp_trialnum = VIP.TrialNum(ii,1:22);
    for m = 1:22
        for n = 1:22
            if temp_trialnum(m)<3 || temp_trialnum(n)<3
                temp_index(m,n) = true;
            else
                temp_index(m,n) = false;
            end
        end
    end
    temp_matrix(temp_index) = nan;
    VIP_Combine.LeverCorr_Reward_Matrix(:,:,ii) = temp_matrix;
    clear temp_matrix temp_trialnum temp_index
end
for jj = 1:size(VIP_2P.LeverCorr_Reward_Matrix,3)
    temp_matrix = VIP_2P.LeverCorr_Reward_Matrix(1:22,1:22,jj);
    temp_trialnum = VIP_2P.TrialNum(jj,1:22);
    for m = 1:22
        for n = 1:22
            if temp_trialnum(m)<3 || temp_trialnum(n)<3
                temp_index(m,n) = true;
            else
                temp_index(m,n) = false;
            end
        end
    end
    temp_matrix(temp_index) = nan;
    VIP_Combine.LeverCorr_Reward_Matrix(:,:,ii+jj) = temp_matrix;
    clear temp_matrix temp_trialnum temp_index
end

Combine.LeverCorr_Reward_Matrix = cat(3,VIP_Combine.LeverCorr_Reward_Matrix,SOM_Combine.LeverCorr_Reward_Matrix);

save('Combine_2PWF_ForBhv','SOM_Combine','VIP_Combine','Combine','-v7.3');

%% Example lever traces
load('Z:\People\Chi\TwoP_IN\VIP\WL_3526642-L\MovAnalysis\WL_3526642-L_CuedMov_All.mat', 'CuedMov_SingleAnimal')
% Use session 1,5,12,20
for ii = [1:21]
    temp_trace{ii} = resample(CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward,1,10);
end
jj = 1;
figure; hold on;
for kk = 1:min(size(temp_trace{jj},2),50)
    subplot(5,10,kk); hold on;
    plot(temp_trace{jj}(10:240,kk),'color',[0.7 0.7 0.7],'linewidth',0.5);
    title(num2str(kk));
    ylim([0.5 2]);
    xlim([1 240]); line([40 40],ylim,'color','k');
end
good_pool = [2,3,4,6,7,8,9,10,11,12,14,16,17,18,19,20,21,23,24,26,27];
ii = 1;
figure; hold on;
for kk = 1:50
    rand_index{ii}(kk,:) = good_pool(randperm(length(good_pool),10));
    selected_trace{ii} = temp_trace{ii}(:,rand_index{ii}(kk,:));
    subplot(5,10,kk); hold on;
    plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
    hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
    ylim([0.5 2]);
    xlim([1 240]); line([40 40],ylim,'color','k');
end
ii = 1; kk = 50;
rand_index_2 = rand_index{ii}(kk,:);
selected_trace{ii} = temp_trace{ii}(:,rand_index_2);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');

jj = 5;
figure; hold on;
for kk = 1:min(size(temp_trace{jj},2),50)
    subplot(5,10,kk); hold on;
    plot(temp_trace{jj}(10:240,kk),'color',[0.7 0.7 0.7],'linewidth',0.5);
    title(num2str(kk));
    ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
end
good_pool = [3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,26,27,28,32,34,35,36,37,39,41,42,43];
ii = 5;
figure; hold on;
for kk = 1:50
    rand_index{ii}(kk,:) = [good_pool(randperm(length(good_pool),8)),24,25];
    selected_trace{ii} = temp_trace{ii}(:,rand_index{ii}(kk,:));
    subplot(5,10,kk); hold on;
    plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
    hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
    ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
end
ii = 5; kk = 40;
rand_index_2 = rand_index{ii}(kk,:);
selected_trace{ii} = temp_trace{ii}(:,rand_index_2);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');

jj = 12;
figure; hold on;
for kk = 1:min(size(temp_trace{jj},2),50)
    subplot(5,10,kk); hold on;
    plot(temp_trace{jj}(10:240,kk),'color',[0.7 0.7 0.7],'linewidth',0.5);
    title(num2str(kk));
    ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
end
good_pool = [2,4,7,8,9,10,11,12,13,14,16];
ii = 12;
figure; hold on;
for kk = 1:50
    rand_index{ii}(kk,:) = good_pool(randperm(length(good_pool),10));
    selected_trace{ii} = temp_trace{ii}(:,rand_index{ii}(kk,:));
    subplot(5,10,kk); hold on;
    plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
    hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
    ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
end
ii = 12; kk = 3;
rand_index_2 = rand_index{ii}(kk,:);
selected_trace{ii} = temp_trace{ii}(:,rand_index_2);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');

jj = 21;
figure; hold on;
for kk = 1:min(size(temp_trace{jj},2),50)
    subplot(5,10,kk); hold on;
    plot(temp_trace{jj}(10:240,kk),'color',[0.7 0.7 0.7],'linewidth',0.5);
    title(num2str(kk));
%     ylim([-1.5 -0.2]);
    xlim([1 240]); line([40 40],ylim,'color','k');
end
good_pool = [3,4,10,18,24,26,27,29,30,31,34,35];
ii = 21;
figure; hold on;
for kk = 1:50
    rand_index{ii}(kk,:) = good_pool(randperm(length(good_pool),10));
    selected_trace{ii} = temp_trace{ii}(:,rand_index{ii}(kk,:));
    subplot(5,10,kk); hold on;
    plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
    hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
%     ylim([-1.5 -0.2]);
    xlim([1 240]); line([40 40],ylim,'color','k');
end
kk = 50;
selected_trace{ii} = temp_trace{ii}(:,rand_index{ii}(kk,:));
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(selected_trace{ii}(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(selected_trace{ii}(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');