%% Plot Figure S8C, cell density, basal forebrain
celldens_ctrl = [70.59073307,46.79696139,81.82644676;106.702035,76.65830896,71.59057815;96.0034281,64.36503523,62.45164616];
celldens_cas3 = [0,0,0;3.554531071,1.478362543,0;0,0,0];
celldens_ctrl = round(celldens_ctrl);
celldens_cas3 = round(celldens_cas3);
figure; set(gcf,'color','w','position',[200 200 200 200]);
markers = {'o','^','s';'o','^','s';'o','^','s'};
hold on;
temp_var_1 = celldens_ctrl;
temp_var_2 = celldens_cas3;
for ii = 1:3
    if ii == 1
        for jj = 1:3
            plot([1,2],[temp_var_1(ii,jj),temp_var_2(ii,jj)],'marker',markers{ii,jj},'color',[0.7,0.7,0.7]);
        end
    elseif ii == 2
        for jj = 1:3
            plot([1,2],[temp_var_1(ii,jj),temp_var_2(ii,jj)],'marker',markers{ii,jj},'color',[0.7,0.7,0.7]);
        end
    elseif ii == 3
        for jj = 1:3
            plot([1,2],[temp_var_1(ii,jj),temp_var_2(ii,jj)],'marker',markers{ii,jj},'color',[0.7,0.7,0.7]);
        end
    end
end
temp_var = [temp_var_1(:),temp_var_2(:)];
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot([1,2],temp_1,'color','k','linewidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color','k','linewidth',2);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'Control','Ablation'})
ylabel('ChAT+ density (cell/mm^2)');
axis square;

%% Plot Figure S8F, cell density, striatum
str_sal = [22.10398798; 13.02462542; 10.63572007; 9.774187162; 8.611613406; 11.97185544; 11.42073663; 10.37764684]; % 5 animals
str_cas = [7.92354758; 12.74566846; 14.41923437; 16.74756007; 13.49899007; 7.896165548; 11.96761341; 6.516861675; 18.47566139]; % 7 animals
figure; hold on; set(gcf,'color','w','position',[200 200 200 300]);
plot(1+((rand(1,length(str_sal(:)))-0.5)*0.25),str_sal(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,length(str_cas(:)))-0.5)*0.25),str_cas(:),'color',[1 0.5 0.5],'marker','o','linestyle','none');
temp_var_11 = nanmean(str_sal,1);
temp_var_12 = nanstd(str_sal,[],1)./sqrt(sum(~isnan(str_sal)));
temp_var_21 = nanmean(str_cas,1);
temp_var_22 = nanstd(str_cas,[],1)./sqrt(sum(~isnan(str_cas)));
bar([1],temp_var_11,'EdgeColor','k','FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor','r','FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color','k','linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color','r','linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'Control','Ablation'});
ylabel('ChAT+ density (cell/mm^2)');
[bootstat_1,~] = bootstrp(20000,@mean,str_sal);
[bootstat_2,~] = bootstrp(20000,@mean,str_cas);
bootstat_diff = bootstat_2-bootstat_1;
p = min(sum(bootstat_diff<0),sum(bootstat_diff>0))/20000;

%% Plot Figure S8G, Casp3 vs Sal
load('FigureS8G.mat');
% Time from cue to rewarded movement onset
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Cue_to_CuedRewardedMov;
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(Cas_index_good),:);
temp_var_2 = temp_var(logical(Sal_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:21
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:21
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 21.5]); xticks([1:2:21]); xlabel('Session');
ylim([0 4]); yticks([0:1:4]); ylabel('Time (s)');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Cue to movement onset');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 3.5;
for ii = 1:length(stage_sessions)
    curr_stage = stage_sessions{ii};
    y = [temp_var_1(:,curr_stage);temp_var_2(:,curr_stage)];
    Drug = [ones(size(temp_var_1(:,curr_stage)));zeros(size(temp_var_2(:,curr_stage)))];
    Sessions = repmat([1:size(y,2)],size(y,1),1);
    Sessions = Sessions(:);
    Animals_test = repmat([1:size(y,1)]',1,size(y,2));
    Animals_test = Animals_test(:);
    Drug = Drug(:);
    y = y(:);
    tbl = table(Drug,Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Drug = nominal(tbl.Drug);
    lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    line([curr_stage(1),curr_stage(end)],[y_pos,y_pos],'color','k')
    text(nanmean(curr_stage),y_pos*1.05,[num2str(pValue(2))],'horizontalalignment','center','fontsize',6);
end

% Time rewarded movement onset to reward
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = CuedRewardedMov_to_Reward;
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(Cas_index_good),:);
temp_var_2 = temp_var(logical(Sal_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:21
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:21
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 21.5]); xticks([1:2:21]); xlabel('Session');
ylim([0.2 0.7]); yticks([0.2:0.1:0.7]); ylabel('Time (s)');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Movement onset to reward');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.65;
for ii = 1:length(stage_sessions)
    curr_stage = stage_sessions{ii};
    y = [temp_var_1(:,curr_stage);temp_var_2(:,curr_stage)];
    Drug = [ones(size(temp_var_1(:,curr_stage)));zeros(size(temp_var_2(:,curr_stage)))];
    Sessions = repmat([1:size(y,2)],size(y,1),1);
    Sessions = Sessions(:);
    Animals_test = repmat([1:size(y,1)]',1,size(y,2));
    Animals_test = Animals_test(:);
    Drug = Drug(:);
    y = y(:);
    tbl = table(Drug,Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Drug = nominal(tbl.Drug);
    lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    line([curr_stage(1),curr_stage(end)],[y_pos,y_pos],'color','k')
    text(nanmean(curr_stage),y_pos*1.05,[num2str(pValue(2))],'horizontalalignment','center','fontsize',6);
end

%% Plot Figure S8I & S8J, the effecs of inactivating axonal activity with eOPN3
load('FigureS8IS8J.mat')
% Figure S8I, heatmap
color_map = cbrewer('div','RdBu',64);
color_map = flipud(color_map);
figure; hold on; set(gcf,'color','w','position',[200 200 200*length(df_f_aligned_CaOn_MeanXTrial_sub) 200])
for curr_session = 1:length(df_f_aligned_CaOn_MeanXTrial_sub)
    subplot(1,length(df_f_aligned_CaOn_MeanXTrial_sub),curr_session);
    temp_matrix = df_f_aligned_CaOn_MeanXTrial_sub{curr_session};
    imagesc(temp_matrix,[-1.5 1.5]); colormap(color_map)
    line([29 29],ylim,'color','k','linestyle',':','linewidth',1);
    xticks([1:28:28*5]);
    xticklabels([-1:1:3]);
    xlabel('Time (s)'); ylabel('ROI (axonal segments)');
end

% Figure S8J, activity trace and level
figure; hold on; set(gcf,'color','w','position',[200 200 200 200])
y_limit = [-0.2 0.6];
for curr_session = 1:length(df_f_aligned_CaOn_MeanXTrial_sub)
    temp = df_f_aligned_CaOn_MeanXTrial_sub{curr_session};
    temp_11 = nanmean(temp,1);
    temp_12 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1)),1));
    h = area([(temp_11-temp_12)',(2*temp_12)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    if curr_session == 2
        set(h(2),'EdgeColor','none','FaceColor','g','FaceAlpha',0.3);
        plot(temp_11,'color','g','linewidth',1);
    else
        set(h(2),'EdgeColor','none','FaceColor',repmat([-0.2+curr_session*0.2],1,3),'FaceAlpha',0.3);
        plot(temp_11,'color',repmat([-0.2+curr_session*0.2],1,3),'linewidth',1);
    end
    xlim([1 length(temp_11)]);
    xticks([1:28:28*5]);
    xticklabels([-1:1:3]);
    xlabel('Time (s)'); ylabel('mean df/f');
    ylim(y_limit);
    line([29 29],ylim,'color','k','linestyle',':');
end
axis square;

figure; hold on; set(gcf,'color','w','position',[200 200 300 200])
temp = df_f_aligned_CaOn_MeanXTrialXTime_sub;
temp_11 = nanmean(temp,1);
temp_12 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1)),1));
yyaxis left;
plot(temp_11,'linewidth',1);
for ii = 1:size(df_f_aligned_CaOn_MeanXTrialXTime_sub,2)
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'linewidth',1);
end
xlim([0.5 3.5]); ylim([0 0.4]);
xticks([1:3]);
xticklabels({'Pre.','2 min p','10 min p','20 min p'});
ylabel('Mean df/f');
temp = df_f_aligned_CaOn_MeanXTrialXTime_sub;
temp = temp./repmat(temp(:,1),1,size(temp,2));
temp_11 = nanmean(temp,1);
temp_12 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1)),1));
yyaxis right;
plot(temp_11,'linewidth',1);
for ii = 1:size(df_f_aligned_CaOn_MeanXTrialXTime_sub,2)
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'linewidth',1);
end
ylabel('Norm. mean df/f'); ylim([0 1])
axis square;
% Stats
aaa = df_f_aligned_CaOn_MeanXTrialXTime_sub;
pvalue(1) = signrank(aaa(:,1),aaa(:,2));
pvalue(2) = signrank(aaa(:,1),aaa(:,3));
pvalue(3) = signrank(aaa(:,2),aaa(:,3));
FDR = mafdr(pvalue,'BHFDR', true);

%% Plot Figure S8K, eOPN3 vs mCherry
load('FigureS8K.mat');
% Time from cue to rewarded movement onset
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Cue_to_CuedRewardedMov;
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(eOPN3_index_good),:);
temp_var_2 = temp_var(logical(mCherry_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:21
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:21
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 21.5]); xticks([1:2:21]); xlabel('Session');
ylim([0 3.2]); yticks([0:1:3]); ylabel('Time (s)');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Cue to movement onset');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 3;
for ii = 1:length(stage_sessions)
    curr_stage = stage_sessions{ii};
    y = [temp_var_1(:,curr_stage);temp_var_2(:,curr_stage)];
    Drug = [ones(size(temp_var_1(:,curr_stage)));zeros(size(temp_var_2(:,curr_stage)))];
    Sessions = repmat([1:size(y,2)],size(y,1),1);
    Sessions = Sessions(:);
    Animals_test = repmat([1:size(y,1)]',1,size(y,2));
    Animals_test = Animals_test(:);
    Drug = Drug(:);
    y = y(:);
    tbl = table(Drug,Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Drug = nominal(tbl.Drug);
    lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    line([curr_stage(1),curr_stage(end)],[y_pos,y_pos],'color','k')
    text(nanmean(curr_stage),y_pos*1.05,[num2str(pValue(2))],'horizontalalignment','center','fontsize',6);
end

% Time rewarded movement onset to reward
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = CuedRewardedMov_to_Reward;
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(eOPN3_index_good),:);
temp_var_2 = temp_var(logical(mCherry_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:21
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:21
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 21.5]); xticks([1:4:21]); xlabel('Session');
ylim([0 1]); yticks([0:0.2:1]); ylabel('Time (s)');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Movement onset to reward');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.8;
for ii = 1:length(stage_sessions)
    curr_stage = stage_sessions{ii};
    y = [temp_var_1(:,curr_stage);temp_var_2(:,curr_stage)];
    Drug = [ones(size(temp_var_1(:,curr_stage)));zeros(size(temp_var_2(:,curr_stage)))];
    Sessions = repmat([1:size(y,2)],size(y,1),1);
    Sessions = Sessions(:);
    Animals_test = repmat([1:size(y,1)]',1,size(y,2));
    Animals_test = Animals_test(:);
    Drug = Drug(:);
    y = y(:);
    tbl = table(Drug,Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Drug = nominal(tbl.Drug);
    lme = fitlme(tbl,'y ~ 1 + Drug + (1|Sessions)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    line([curr_stage(1),curr_stage(end)],[y_pos,y_pos],'color','k')
    text(nanmean(curr_stage),y_pos*1.05,[num2str(pValue(2))],'horizontalalignment','center','fontsize',6);
end

