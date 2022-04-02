%% Plot Figure 1B, correct rate
load('Figure1B1D1E.mat')
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.CR;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:22
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 22.5]);xticks([5:5:20]);xlabel('Session');
ylim([0 1]); yticks([0:0.2:1]); yticklabels([0:20:100]); ylabel('% correct');
set(gca,'XColor','k','YColor','k');
line([2.5 2.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([8.5 8.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([16.5 16.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,0.3,['Session: ' num2str(pValue(2))],'fontsize',6);

%% Plot Figure 1C
% Example: WL_3526642-L, session 1, 5, 12, 21
load('Figure1C.mat');
% session 1
curr_random_index = [2,3,4,7,11,16,18,21,23,26];
example_session_1 = LeverTraces{1}(:,curr_random_index);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(example_session_1(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(example_session_1(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
% session 5
curr_random_index = [2,7,8,11,12,13,14,16,24,25];
example_session_5 = LeverTraces{5}(:,curr_random_index);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(example_session_5(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(example_session_5(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
% session 12
curr_random_index = [2,4,7,8,9,11,12,13,14,16];
example_session_12 = LeverTraces{12}(:,curr_random_index);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(example_session_12(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(example_session_12(10:240,:),2),'color','r','linewidth',1);
ylim([0.5 2]);xlim([1 240]); line([40 40],ylim,'color','k');
% session 21
curr_random_index = [3,4,18,24,26,27,30,31,34,35];
example_session_21 = LeverTraces{21}(:,curr_random_index);
figure; set(gcf,'color','w','position',[200,200,300,200]); hold on;
plot(example_session_21(10:240,:),'color',[0.7 0.7 0.7],'linewidth',0.5);
hold on; plot(nanmean(example_session_21(10:240,:),2),'color','r','linewidth',1);
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

%% Plot Figure 1D, correlation matrix
load('Figure1B1D1E.mat');
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = nanmean(Combine.LeverCorr_Reward_Matrix,3); % already excluded session with <3 trials
imagesc(temp_var,[0 0.3]); colormap hot;
xlim([0.5 22.5]); ylim([0.5 22.5]);
set(gca,'YDir','reverse');
xlabel('Session'); ylabel('Session');
axis square

%% Plot Figure 1E, lever correlation
load('Figure1B1D1E.mat')
% Within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.LeverCorr_Reward_within;
temp_var(Combine.TrialNum<3) = nan;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:22
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 22.5]);xticks([5:5:20]);xlabel('Session');
ylim([0 0.35]); yticks([0:0.1:0.3]); ylabel('Lever correlation');
set(gca,'XColor','k','YColor','k');
line([2.5 2.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([8.5 8.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([16.5 16.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
title('Within sessions');
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,0.3,['Session: ' num2str(pValue(2))],'fontsize',6);

% Across sessions
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.LeverCorr_Reward_across;
trial_index = (Combine.TrialNum(:,1:end-1)<3)+(Combine.TrialNum(:,2:end)<3);
temp_var(logical(trial_index)) = nan;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:21
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 21.5]);xticks([5,13,20]); xticklabels({'5-6','13-14','20-21'}); xlabel('Session');
ylim([0 0.3]); yticks([0:0.1:0.3]); ylabel('Lever correlation');
set(gca,'XColor','k','YColor','k');
line([1.5 1.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([7.5 7.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([15.5 15.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
axis square
title('Across sessions');
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,0.3,['Session: ' num2str(pValue(2))],'fontsize',6);

