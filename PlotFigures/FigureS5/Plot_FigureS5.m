%% Plot Figure S5B, VIP-IN activty suppressed by CNO and hM4Di, trace
load('FigureS5BS5C.mat','MvmOnsetAligned_df_all','-mat');
figure; set(gcf,'color','w','position',[100 -100 1200 200]); hold on;
for ii = 1:4
    subplot(1,4,ii);hold on;
    temp = MvmOnsetAligned_df_all{ii};
    temp_1 = nanmean(temp,1);
    temp_2 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor','k','FaceAlpha',0.2);
    plot(temp_1,'k','linewidth',1);
    ylim([-0.12 0.3]);
    line([8 8],ylim,'color','k','linestyle',':');
    xlim([1 36]); xticks([8 22 36]); xticklabels({'0','1','2'});
    xlabel('Time (s)'); ylabel('Mean df/f');
    title(sessions_order{ii});
    axis square;
end

%% Plot Figure S5C, VIP-IN activty suppressed by CNO and hM4Di, activity level
figure; set(gcf,'color','w','position',[100 -100 300 200]); hold on;
for ii = 1:size(Mean_PostMvmOnset,1)
    plot(Mean_PostMvmOnset(ii,:),'marker','o','markersize',6,'color',[0.7 0.7 0.7]);
end
plot(nanmean(Mean_PostMvmOnset),'marker','.','markersize',24,'color','k','linewidth',2);
temp_1 = nanmean(Mean_PostMvmOnset);
temp_2 = nanstd(Mean_PostMvmOnset,[],1)./sqrt(sum(~isnan(nanmean(Mean_PostMvmOnset(:,1)))));
for ii = 1:4
    line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color','k','linewidth',2);
end
xlim([0.5 4.5]); xticklabels({'Before','30 min','60 min','120 min'});
ylabel('Mean df/f');
% Stats
for ii = 2:4
    temp_boot = Mean_PostMvmOnset(:,ii)-Mean_PostMvmOnset(:,1);
    [bootstat,~] = bootstrp(20000,@mean,temp_boot);
    if prctile(bootstat,100-0.001*50)<0
        text(ii,0.6,'***','HorizontalAlignment','center');
    elseif prctile(bootstat,100-0.01*50)<0
        text(ii,0.6,'**','HorizontalAlignment','center');
    elseif prctile(bootstat,100-0.05*50)<0
        text(ii,0.6,'*','HorizontalAlignment','center');
    end
end

%% Plot Figure S5E, effectve window of repeated injection of CNO
load('FigureS5E.mat');
MovOnset_Frame = 8;
Baseline_Frame = [3:5];
% Example activity traces
g_n = 10; % Animal 2, neuron 10 in FOV 1, GC6f only
figure; hold on; set(gcf,'color','w','pos',[50 50 400 150]);
hold on;
for ii = 1:3
    subplot(1,3,ii); hold on;
    temp_trace = df_f_MovOnset_Alinged{1,ii}{g_n};
    temp_trace = temp_trace-repmat(nanmean(temp_trace(Baseline_Frame,:)),36,1);
    temp_trace = temp_trace';
    temp_1 = nanmean(temp_trace);
    temp_2 = nanstd(temp_trace)/sqrt(size(temp_trace,1));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[90,110,60]/255,'FaceAlpha',0.2);
    plot(temp_1,'color',[90,110,60]/255,'linewidth',1);
    ylim([-1 4.5]); ylabel('Mean df/f');
    xlim([0 36]); xticks([8,22,36]); xticklabels([0,1,2]); xlabel('Time (s)');
    line([8 8],ylim,'color','k','linestyle',':');
    axis square;
end
g_n = 27; % Animal 2, neuron 27 in FOV 1, GC6f only
figure; hold on; set(gcf,'color','w','pos',[50 50 400 150]);
hold on;
for ii = 1:3
    subplot(1,3,ii); hold on;
    temp_trace = df_f_MovOnset_Alinged{1,ii}{g_n};
    temp_trace = temp_trace-repmat(nanmean(temp_trace(Baseline_Frame,:)),36,1);
    temp_trace = temp_trace';
    temp_1 = nanmean(temp_trace);
    temp_2 = nanstd(temp_trace)/sqrt(size(temp_trace,1));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[90,110,60]/255,'FaceAlpha',0.2);
    plot(temp_1,'color',[90,110,60]/255,'linewidth',1);
    ylim([-1 4.5]); ylabel('Mean df/f');
    xlim([0 36]); xticks([8,22,36]); xticklabels([0,1,2]); xlabel('Time (s)');
    line([8 8],ylim,'color','k','linestyle',':');
    axis square;
end
g_n = 16; % Animal 2, neuron 10 in FOV 1, co-expression
figure; hold on; set(gcf,'color','w','pos',[50 50 400 150]);
hold on;
for ii = 1:3
    subplot(1,3,ii); hold on;
    temp_trace = df_f_MovOnset_Alinged{1,ii}{g_n};
    temp_trace = temp_trace-repmat(nanmean(temp_trace(Baseline_Frame,:)),36,1);
    temp_trace = temp_trace';
    temp_1 = nanmean(temp_trace);
    temp_2 = nanstd(temp_trace)/sqrt(size(temp_trace,1));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[210,110,45]/255,'FaceAlpha',0.2);
    plot(temp_1,'color',[210,110,45]/255,'linewidth',1);
    ylim([-1 4.5]); ylabel('Mean df/f');
    xlim([0 36]); xticks([8,22,36]); xticklabels([0,1,2]); xlabel('Time (s)');
    line([8 8],ylim,'color','k','linestyle',':');
    axis square;
end    
g_n = 3; % Animal, neuron 3 in FOV 2, co-expression
figure; hold on; set(gcf,'color','w','pos',[50 50 400 150]);
hold on;
for ii = 1:3
    subplot(1,3,ii); hold on;
    temp_trace = df_f_MovOnset_Alinged{2,ii}{g_n};
    temp_trace = temp_trace-repmat(nanmean(temp_trace(Baseline_Frame,:)),36,1);
    temp_trace = temp_trace';
    temp_1 = nanmean(temp_trace);
    temp_2 = nanstd(temp_trace)/sqrt(size(temp_trace,1));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[210,110,45]/255,'FaceAlpha',0.2);
    plot(temp_1,'color',[210,110,45]/255,'linewidth',1);
    ylim([-1 4.5]); ylabel('Mean df/f');
    xlim([0 36]); xticks([8,22,36]); xticklabels([0,1,2]); xlabel('Time (s)');
    line([8 8],ylim,'color','k','linestyle',':');
    axis square;
end    

% Activity level
figure; hold on; set(gcf,'color','w','pos',[200,200,200,200]);
temp = cell2mat(Green_Post_conc_3); % GC6f only
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[90,110,60]/255,'FaceAlpha',0.2);
plot(temp_1,'marker','o','color',[90,110,60]/255,'linewidth',1);
temp = cell2mat(Orange_Post_conc_3); % co-expression of GC6f and hM4Di
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[210,110,45]/255,'FaceAlpha',0.2);
plot(temp_1,'marker','o','color',[210,110,45]/255,'linewidth',1);
xlim([0.5 3.5]); xlabel('Session'); ylabel('Mean df/f'); axis square;
xticks([1:1:3]);
% Stats
temp_1 = cell2mat(Green_Post_conc_3);
temp_2 = cell2mat(Orange_Post_conc_3);
[p,h] = ranksum(temp_1(:,1),temp_2(:,1))
[p,h] = ranksum(temp_1(:,2),temp_2(:,2))
[p,h] = ranksum(temp_1(:,3),temp_2(:,3))

%% Plot Figure S5F, activity distribution, randomly selected 50% animals every sweep
load('FigureS5F.mat');
edges = [0.09:0.01:0.32];
figure; set(gcf,'position',[200,200,200,200]); hold on;
histogram(combo_matrix(:,2),edges,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.3,'Normalization','probability');
histogram(combo_matrix(:,1),edges,'edgecolor','none','facecolor',color_value,'facealpha',0.3,'Normalization','probability');
xlim([0.08 0.33]); xlabel('df/f'); ylabel('Prob.')
ylim([0 0.2]); axis square;
[h,p] = kstest2(combo_matrix(:,1),combo_matrix(:,2));
text(0.2,0.17,[num2str(p)],'horizontalalignment','center','fontsize',8);

%% Plot Figure S5G, basic behavior measurement
load('FigureS5G.mat');
% fraction of responded trials
temp_var_1 = mCherry_rspdfraction;
temp_var_2 = hM4Di_rspdfraction;
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
ylabel('Fraction (%)'); title('Responded trials'); axis square;
% Stats
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
text(1.5,1,['hM4Di: ' num2str(pValue(2))],'fontsize',6);

% Time spent in lever-press movements
temp_var_1 = mCherry_mvmfraction;
temp_var_2 = hM4Di_mvmfraction;
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
ylim([0 0.6]); yticks([0:0.1:0.6]); yticklabels([0:10:60]);
ylabel('Fraction (%)'); title('Time in lever-press movements'); axis square;
% Stats
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
text(1.5,0.55,['hM4Di: ' num2str(pValue(2))],'fontsize',6);

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
ylim([0 15]); ylabel('Speed (mm/s)'); title('Movement speed');
axis square;
% Stats
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
