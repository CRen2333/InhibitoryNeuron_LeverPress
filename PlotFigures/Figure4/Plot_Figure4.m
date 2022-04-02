%% VIP
clear;
load('Figure4.mat','VIP','-mat');
Animals = VIP.Animals;
color_value = VIP.color_value;
color_line = VIP.color_line;
All_cell_stage_conc = VIP.All_cell_stage_conc;
Delta_df_against_N = VIP.Delta_df_against_N;
all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages = VIP.all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages;
mycmap = VIP.mycmap;
MovOnset_Frame = VIP.MovOnset_Frame;
All_Delta_df_against_N_CellTag = VIP.All_Delta_df_against_N_CellTag;

%% Plot Figure 4B, activty heatmap
figure; hold on; set(gcf,'color','w','position',[1500 50 800 500]);
Stages = {'Naive','Early','Middle','Late'};
for ii = 1:4
    subplot(1,4,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-0.7 2.6]); colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Stages{ii});
end

%% Plot Figure 4C, average activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
ylim([-0.1 0.6]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
ylabel('Mean df/f'); axis square;

%% Plot Figure 4D, averaged activity level
stages = {'Naive','Early','Middle','Late'};
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = All_cell_stage_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,0.35]); xticks([1:4]); xticklabels(stages); ylabel('Mean df/f');
yticks([0:0.1:0.3]);
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.25,'p < 0.0001','color','k','HorizontalAlignment','center')
axis square

% Compare to naive stage
tbl.Sessions = nominal(tbl.Sessions);
glme = fitglme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[gbeta,gbetanames,gstats] = fixedEffects(glme);
PValues = gstats.pValue(2:end);
[FDR] = mafdr(PValues,'BHFDR', true);

%% Plot Figure 4E, changes of activity compare to naive stage, distribution
figure; set(gcf,'color','w','position',[100,100,600,200]); hold on;
edges = [-2:0.1:2];
for ii = 2:4
    subplot(1,4,ii);
    temp_n = histcounts(Delta_df_against_N(:,ii),edges,'Normalization','probability');
    histogram(Delta_df_against_N(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(ii,:));
    ylim([0 0.25]); xlim([edges(1) edges(end)]);
    line([0 0],ylim,'color','k','linestyle',':');
    xlabel('delta df/f'); ylabel('Prob.'); axis square;
end

%% Plot Figure 4F, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
pie([sum(All_Delta_df_against_N_CellTag(:,4)==1),sum(All_Delta_df_against_N_CellTag(:,4)==-1),sum(All_Delta_df_against_N_CellTag(:,4)==0)]);
axis square;

%% SOM
clear;
load('Figure4.mat','SOM','-mat');
Animals = SOM.Animals;
color_value = SOM.color_value;
color_line = SOM.color_line;
All_cell_stage_conc = SOM.All_cell_stage_conc;
Delta_df_against_N = SOM.Delta_df_against_N;
all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages = SOM.all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages;
mycmap = SOM.mycmap;
MovOnset_Frame = SOM.MovOnset_Frame;
All_Delta_df_against_N_CellTag = SOM.All_Delta_df_against_N_CellTag;

%% Plot Figure 4G, activty heatmap
figure; hold on; set(gcf,'color','w','position',[50 50 800 500]);
Stages = {'Naive','Early','Middle','Late'};
for ii = 1:4
    subplot(1,4,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
    imagesc(temp,[-0.8 1.9]); colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Stages{ii});
end

%% Plot Figure 4H, average activity traces
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
ylim([-0.03 0.2]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
yticks([0:0.1:0.2]);
ylabel('Mean df/f'); axis square

%% Plot Figure 4I, averaged activity level
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = All_cell_stage_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,0.15]); xticks([1:4]); xticklabels(stages); ylabel('Mean df/f');
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.05,'p = 0.0007','color','k','HorizontalAlignment','center')
axis square

% Compare to naive stage
tbl.Sessions = nominal(tbl.Sessions);
glme = fitglme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[gbeta,gbetanames,gstats] = fixedEffects(glme);
PValues = gstats.pValue(2:end);
[FDR] = mafdr(PValues,'BHFDR', true);

%% Plot Figure 4J, changes of activity compare to naive stage, distribution
figure; set(gcf,'color','w','position',[100,100,600,200]); hold on;
edges = [-1.2:0.1:1.2];
for ii = 2:4
    subplot(1,4,ii);
    temp_n = histcounts(Delta_df_against_N(:,ii),edges,'Normalization','probability');
    histogram(Delta_df_against_N(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(ii,:));
    ylim([0 0.4]);
    xlim([edges(1) edges(end)]);
    line([0 0],ylim,'color','k','linestyle',':');
    xlabel('delta df/f'); ylabel('Prob.'); axis square;    
end

%% Plot Figure 4K, significantly changed fraction
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
pie([sum(All_Delta_df_against_N_CellTag(:,ii)==1),sum(All_Delta_df_against_N_CellTag(:,ii)==-1),sum(All_Delta_df_against_N_CellTag(:,ii)==0)]);
axis square;

