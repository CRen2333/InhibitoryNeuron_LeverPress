%% Plot Figure 5, Activity (5B-5E)
load('Figure5_Activity.mat');

%% Plot Figure 5B, Activity heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 400 300]);
subplot(1,2,1);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_m);
[~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
temp = temp(I,:);
imagesc(temp,[-0.8 2]);colormap(mycmap);
line([8,8],ylim,'color','k','linestyle',':');
ylabel('# of cells');
xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
title('mCherry');
subplot(1,2,2);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_h);
[~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
temp = temp(I,:);
imagesc(temp,[-0.8 2]);colormap(mycmap);
line([8,8],ylim,'color','k','linestyle',':');
ylabel('# of cells');
xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
title('hM4Di');

%% Plot Figure 5C, Average activity traces during movement
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_m);
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5 0.5 0.5],'linewidth',1);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_h);
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value,'FaceAlpha',0.3);
plot(temp_1,'color',color_value,'linewidth',2);
ylim([-0.05 0.45]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square

%% Plot Figure 5D, Activity level during movements
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = All_cell_stage_conc_m;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
bar([1],temp_1,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_1(1)-temp_2(1),temp_1(1)+temp_2(1)],'color',[0.5 0.5 0.5],'LineWidth',1);
temp_var = All_cell_stage_conc_h;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
bar([2],temp_1,'EdgeColor',color_value,'FaceColor','none','linewidth',1,'barwidth',0.7);
line([2,2],[temp_1(1)-temp_2(1),temp_1(1)+temp_2(1)],'color',color_value,'LineWidth',1);
xlim([0.4,2.6]); ylim([0,0.3]); xticks([1:2]); xticklabels({'mCherry','hM4Di'}); ylabel('Mean df/f');
axis square;
% Stats
for ii = 1:length(Animals)
    NeuroNum{ii,1} = length(cell2mat(Post_MovOnset_Aligend_df_field_stage{ii}));
    Animals_tag{ii,1} = repmat(ii,NeuroNum{ii,1},1);
    hM4Di_test{ii,1} = repmat(hM4Di_index(ii),NeuroNum{ii,1},1);
end
Drug = cell2mat(hM4Di_test);
Animals_test = cell2mat(Animals_tag);
Animals_test = Animals_test(:);
y = cell2mat(Post_MovOnset_Aligend_df_stage');
tbl = table(Drug,Animals_test,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;    
line([0.9,2.1],[0.25 0.25],'color','k')
text(1.5,0.27,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

%% Plot Figure 5E, Activity distribution during movements
edges = [-1.4:0.1:2];
figure; set(gcf,'position',[200,200,200,200]); hold on;
histogram(All_cell_stage_conc_m,edges,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.3,'Normalization','probability');
histogram(All_cell_stage_conc_h,edges,'edgecolor','none','facecolor',color_value,'facealpha',0.3,'Normalization','probability');
xlim([-1.4 2]); xlabel('Mean df/f'); ylabel('Prob.')
ylim([0 0.36]); axis square;

%% Plot Figure 5F, Behavior, naive stage (5F)
load('Figure5_Behavior.mat');
% Fraction of rewarded trials
temp_var = VIPSOM_2P.CR;
temp_var_1 = temp_var(mCherry_index,1:2);
temp_var_2 = temp_var(hM4Di_index,1:2);
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
ylim([0 1]); yticks([0:0.2:1]); yticklabels([0:20:100]);
ylabel('% correct'); title('Rewarded trials'); axis square;
% Stats
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
line([1,2],[0.95 0.95],'color','k')
text(1.5,0.97,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);

% lever correlation, within sessions
temp_var = VIPSOM_2P.LeverCorr_Reward_within;
temp_var(VIPSOM_2P.TrialNum<3) = nan;
temp_var_1 = temp_var(mCherry_index,1:2);
temp_var_2 = temp_var(hM4Di_index,1:2);
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
ylim([-0.22 0.4]); yticks([-0.2:0.2:0.4]);
ylabel('Lever correlation'); title('Winthin sessions'); axis square;
% Stats
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
line([1,2],[0.37 0.37],'color','k')
text(1.5,0.38,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);

% Cue to movement onset
temp_var = VIPSOM_2P.C_CRM;
temp_var(VIPSOM_2P.TrialNum<3) = nan;
temp_var_1 = temp_var(mCherry_index,1:2);
temp_var_2 = temp_var(hM4Di_index,1:2);
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
ylim([0 5]); yticks([0:1:5]);
ylabel('Time (s)'); title('Cue to movement onset'); axis square;
% Stats
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
line([1,2],[4.5 4.5],'color','k')
text(1.5,4.7,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);

% Movement onset to reward
temp_var = VIPSOM_2P.RwdMVM_Rwd;
temp_var(VIPSOM_2P.TrialNum<3) = nan;
temp_var_1 = temp_var(mCherry_index,1:2);
temp_var_2 = temp_var(hM4Di_index,1:2);
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
ylim([0 2.5]); yticks([0:0.5:2.5]);
ylabel('Time (s)'); title('Movement onset to reward'); axis square;
% Stats
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
line([1,2],[2.2 2.2],'color','k')
text(1.5,2.3,['hM4Di: ' num2str(pValue(2))],'horizontalalignment','center','fontsize',6);

%% Figure 5G, learning effect on behavior performance
load('Figure5G.mat');
% rewarded trials in first 10 trials
temp_var = CR_first10;
temp_var_1 = temp_var(mCherry_index,:);
temp_var_2 = temp_var(hM4Di_index,:);
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
ylim([0 1]); yticks([0:0.2:1]); yticklabels([0:20:100]);
ylabel('% correct'); title('Rewarded trials'); axis square;
% Stats
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
line([0.9,2.1],[0.9 0.9],'color','k')
text(1.5,0.95,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

% Time from movement onset to reward
temp_var = RwdMVM_Rwd_first10;
temp_var_1 = temp_var(mCherry_index,:);
temp_var_2 = temp_var(hM4Di_index,:);
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
ylim([0 3]); yticks([0:1:3]); ylabel('Time (s)');
title('Movement onset to reward'); axis square;
% Stats
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
line([0.9,2.1],[2.75 2.75],'color','k')
text(1.5,3,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);
