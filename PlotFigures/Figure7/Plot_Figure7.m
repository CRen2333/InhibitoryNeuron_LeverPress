%% VIP, naive stage + AChR antagonists (n&m)
load('Figure7_VIP_Nai_Ant_Cocktail.mat');

%% Plot Figure 7B, activity heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
for ii = 1:2
    subplot(1,2,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-1 2.5]);
    colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end

%% Plot Figure 7C, averaged activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 1;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
ylim([-0.08 0.38]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square;

%% Plot Figure 7D, activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
temp_var = All_cell_Pharm_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
area([1.5,2.5],[1,1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
area([1.5,2.5],[-1,-1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
plot(temp_1,'color',color_value(1,:),'LineWidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_value(1,:),'LineWidth',2);
end
xlim([0.5,2.5]); xticks([1:2]); xticklabels(Pharms); ylabel('Mean df/f'); axis square
ylim([0,0.18]);
% random effect modal
Task = repmat([0:1],size(temp_var,1),1);
Task = Task(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
Neurons = repmat([1:size(temp_var,1)],1,2);
Neurons = Neurons(:);
y = temp_var(:);
tbl = table(Animals_test,Task,y,Neurons);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Neurons = nominal(tbl.Neurons);
tbl.Task = nominal(tbl.Task);
% lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)');
lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)+(Task-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,0.15,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square

%% Plot Figure 7E, distribution of activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
edges = [-2:0.1:2];
temp_n = histcounts(Delta_df_against_C(:,2),edges,'Normalization','probability');
histogram(Delta_df_against_C(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(1,:));
ylim([0 0.3]); xlim([edges(1) edges(end)]);
line([0 0],ylim,'color','k','linestyle',':');
xlabel('delta df/f'); ylabel('Prob.'); axis square;

%% Plot Figure 7F, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
pie([sum(All_Delta_df_against_C_CellTag(:,ii)==1),sum(All_Delta_df_against_C_CellTag(:,ii)==-1),sum(All_Delta_df_against_C_CellTag(:,ii)==0)]);
axis square; axis off;

%% SOM, naive stage + AChR antagonists (n&m)
load('Figure7_SOM_Nai_Ant_Cocktail.mat');

%% Plot Figure 7G, activity heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
for ii = 1:2
    subplot(1,2,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-0.8 2.5]);
    colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end

%% Plot Figure 7H, averaged activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 1;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
ylim([-0.05 0.3]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square;

%% Plot Figure 7I, activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
temp_var = All_cell_Pharm_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
area([1.5,2.5],[1,1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
area([1.5,2.5],[-1,-1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
plot(temp_1,'color',color_value(1,:),'LineWidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_value(1,:),'LineWidth',2);
end
xlim([0.5,2.5]); xticks([1:2]); xticklabels(Pharms); ylabel('Mean df/f'); axis square
ylim([0,0.15]);
% random effect modal
Task = repmat([0:1],size(temp_var,1),1);
Task = Task(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
Neurons = repmat([1:size(temp_var,1)],1,2);
Neurons = Neurons(:);
y = temp_var(:);
tbl = table(Animals_test,Task,y,Neurons);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Neurons = nominal(tbl.Neurons);
tbl.Task = nominal(tbl.Task);
% lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)');
lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)+(Task-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,0.15,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square

%% Plot Figure 7J, distribution of activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
edges = [-2:0.1:2];
temp_n = histcounts(Delta_df_against_C(:,2),edges,'Normalization','probability');
histogram(Delta_df_against_C(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(1,:));
ylim([0 0.35]); xlim([edges(1) edges(end)]);
line([0 0],ylim,'color','k','linestyle',':');
xlabel('delta df/f'); ylabel('Prob.'); axis square;

%% Plot Figure 7K, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
pie([sum(All_Delta_df_against_C_CellTag(:,ii)==1),sum(All_Delta_df_against_C_CellTag(:,ii)==-1),sum(All_Delta_df_against_C_CellTag(:,ii)==0)]);
axis square; axis off;

%% VIP, naive stage + nAChR antagonist
load('Figure7_VIP_Nai_nAnt_Cocktail.mat');

%% Plot Figure 7L, activity heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
for ii = 1:2
    subplot(1,2,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-1 2]);
    colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end

%% Plot Figure 7M, averaged activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 1;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
ylim([-0.03 0.3]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square;

%% Plot Figure 7N, activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
temp_var = All_cell_Pharm_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
area([1.5,2.5],[1,1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
area([1.5,2.5],[-1,-1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
plot(temp_1,'color',color_value(1,:),'LineWidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_value(1,:),'LineWidth',2);
end
xlim([0.5,2.5]); xticks([1:2]); xticklabels(Pharms); ylabel('Mean df/f'); axis square
ylim([0,0.2]);
% random effect modal
Task = repmat([0:1],size(temp_var,1),1);
Task = Task(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
Neurons = repmat([1:size(temp_var,1)],1,2);
Neurons = Neurons(:);
y = temp_var(:);
tbl = table(Animals_test,Task,y,Neurons);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Neurons = nominal(tbl.Neurons);
tbl.Task = nominal(tbl.Task);
% lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)');
lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)+(Task-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,0.15,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square

%% Plot Figure 7O, distribution of activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
edges = [-2:0.1:2];
temp_n = histcounts(Delta_df_against_C(:,2),edges,'Normalization','probability');
histogram(Delta_df_against_C(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(1,:));
ylim([0 0.25]); xlim([edges(1) edges(end)]);
line([0 0],ylim,'color','k','linestyle',':');
xlabel('delta df/f'); ylabel('Prob.'); axis square;

%% Plot Figure 7P, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
pie([sum(All_Delta_df_against_C_CellTag(:,ii)==1),sum(All_Delta_df_against_C_CellTag(:,ii)==-1),sum(All_Delta_df_against_C_CellTag(:,ii)==0)]);
axis square; axis off;

%% VIP, naive stage + mAChR antagonist
load('Figure7_VIP_Nai_mAnt_Cocktail.mat');

%% Plot Figure 7Q, activity heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
for ii = 1:2
    subplot(1,2,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-1 2]);
    colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end

%% Plot Figure 7R, averaged activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 1;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
ylim([-0.05 0.45]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square;

%% Plot Figure 7S, activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
temp_var = All_cell_Pharm_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
area([1.5,2.5],[1,1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
area([1.5,2.5],[-1,-1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
plot(temp_1,'color',color_value(1,:),'LineWidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_value(1,:),'LineWidth',2);
end
xlim([0.5,2.5]); xticks([1:2]); xticklabels(Pharms); ylabel('Mean df/f'); axis square
ylim([0,0.35]);
% random effect modal
Task = repmat([0:1],size(temp_var,1),1);
Task = Task(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
Neurons = repmat([1:size(temp_var,1)],1,2);
Neurons = Neurons(:);
y = temp_var(:);
tbl = table(Animals_test,Task,y,Neurons);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Neurons = nominal(tbl.Neurons);
tbl.Task = nominal(tbl.Task);
% lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)');
lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)+(Task-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,0.15,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square

%% Plot Figure 7T, distribution of activity level
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
edges = [-2:0.1:2];
temp_n = histcounts(Delta_df_against_C(:,2),edges,'Normalization','probability');
histogram(Delta_df_against_C(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(1,:));
ylim([0 0.2]); xlim([edges(1) edges(end)]);
line([0 0],ylim,'color','k','linestyle',':');
xlabel('delta df/f'); ylabel('Prob.'); axis square;

%% Plot Figure 7U, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
pie([sum(All_Delta_df_against_C_CellTag(:,ii)==1),sum(All_Delta_df_against_C_CellTag(:,ii)==-1),sum(All_Delta_df_against_C_CellTag(:,ii)==0)]);
axis square; axis off;


