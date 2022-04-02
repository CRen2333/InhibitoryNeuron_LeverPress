%% Plot Figure 6C, cholinergic activity
% M1
load('Figure6C_M1.mat');
color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];   
color_line = nanmean(color_value(2:3,:));
stages = {'Naive','Early','Middle','Late'};
% Activity traces
clear temp_mean
for curr_animal = 1:length(Animals)
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,1};
    temp_mean{1}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,2};
    temp_mean{2}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,3};
    temp_mean{3}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,4};
    temp_mean{4}(curr_animal,:) = nanmean(temp,1);
end
% across animals
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = temp_mean{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
xlim([1 76]); ylim([-0.1 0.7])
xticks([16,46,76]);xticklabels({'0','1','2'});
xlabel('Time (s)');
ylabel('Mean df/f');
line([16 16], ylim, 'color','k','linestyle',':','linewidth',1);
axis square;
% Activity level
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = Post_MovOnset_Aligend_df_stage_sub_Mean_2;
for curr_animal = 1:size(temp_var,1)
    plot(temp_var(curr_animal,:),'color',[0.5 0.5 0.5],'LineWidth',0.5);
end
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,1]); xticklabels(stages); ylabel('Mean df/f');
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,4);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.8,'p < 0.01','color','k','fontsize',8,'HorizontalAlignment','center');
xlabel('Stages'); ylabel('Mean df/f');
axis square

% S1 and PPC
load('Figure6C_S1HL_PPC.mat')
color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];   
color_line = nanmean(color_value(2:3,:));
stages = {'Naive','Early','Middle','Late'};
% S1HL
% Activity traces
clear temp_mean
curr_field = 1;
for curr_animal = 1:length(Animals)
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,1};
    temp_mean{1}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,2};
    temp_mean{2}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,3};
    temp_mean{3}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,4};
    temp_mean{4}(curr_animal,:) = nanmean(temp,1);
end
% across animals
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = temp_mean{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
xlim([1 76]); ylim([-0.2 1.2]); yticks([0:0.4:1.2]);
xticks([16,46,76]);xticklabels({'0','1','2'});
xlabel('Time (s)');
ylabel('Mean df/f');
line([16 16], ylim, 'color','k','linestyle',':','linewidth',1);
axis square;
% Activity level
curr_field = 1;
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{curr_field};
for curr_animal = 1:size(temp_var,1)
    plot(temp_var(curr_animal,:),'color',[0.5 0.5 0.5],'LineWidth',0.5);
end
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,2]); xticklabels(stages); ylabel('Mean df/f');
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,4);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,1.5,['p = ' num2str(pValue(2), '%.4f')],'color','k','fontsize',8,'HorizontalAlignment','center')
axis square
title('S1HL','fontsize',10);
%PPC
% Activity traces
clear temp_mean
curr_field = 2;
for curr_animal = 1:length(Animals)
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,1};
    temp_mean{1}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,2};
    temp_mean{2}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,3};
    temp_mean{3}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,4};
    temp_mean{4}(curr_animal,:) = nanmean(temp,1);
end
% across animals
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = temp_mean{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
xlim([1 76]); ylim([-0.1 0.65]); yticks([0:0.2:0.6]);
xticks([16,46,76]);xticklabels({'0','1','2'});
xlabel('Time (s)');
ylabel('Mean df/f');
line([16 16], ylim, 'color','k','linestyle',':','linewidth',1);
axis square;
% Activity level
curr_field = 2;
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{curr_field};
for curr_animal = 1:size(temp_var,1)
    plot(temp_var(curr_animal,:),'color',[0.5 0.5 0.5],'LineWidth',0.5);
end
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([-0.1 0.7]); xticklabels(stages); ylabel('Mean df/f');
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,4);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.5,['p = ' num2str(pValue(2), '%.4f')],'color','k','fontsize',8,'HorizontalAlignment','center')
axis square
title('PPC','fontsize',10);

%% Plot Figure 6E, correct rate, Casp3 vs Sal
load('Figure6E6F.mat');

figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = CR;
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
xlim([0.5 21.5]); xticks([1:4:21]); xlabel('Session');
yticks([0:0.2:1]); yticklabels([0:20:100]); ylabel('% correct');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Rewarded trials');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.2;
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

%% Plot Figure 6F, lever correlation, Casp3 vs Sal
load('Figure6E6F.mat');
% within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Corr_Reward_within;
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
xlim([0.5 21.5]); xticks([1:4:21]); xlabel('Session');
ylim([0 0.45]); yticks([0:0.1:0.4]); ylabel('Lever correlation');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Within sessions');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.4;
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
% across sessions
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Corr_Reward_across;
temp_var_1 = temp_var(logical(Cas_index_good),:);
temp_var_2 = temp_var(logical(Sal_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:20
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:20
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 20.5]); xticks([1,7,13,19]); xticklabels({'1-2','7-8','13-14','19-20'}); xlabel('Session');
ylim([0 0.4]); yticks([0:0.1:0.4]); ylabel('Lever correlation');
line([1.5 1.5],ylim,'color','k','linestyle',':');
line([7.5 7.5],ylim,'color','k','linestyle',':');
line([14.5 14.5],ylim,'color','k','linestyle',':');
title('Across sessions');
% Stats
stage_sessions = {[1],[2:7],[8:14],[15:20]};
y_pos = 0.35;
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

%% Plot Figure 6H, correct rate, eOPN3 vs mCherry
load('Figure6H6I.mat');
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = CR;
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
yticks([0:0.2:1]); yticklabels([0:20:100]); ylabel('% correct');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Rewarded trials');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.2;
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

%% Plot Figure 6I, lever correlation, eOPN3 vs mCherry
load('Figure6H6I.mat');
% within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Corr_Reward_within;
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
ylim([0 0.4]); yticks([0:0.1:0.4]); ylabel('Lever correlation');
line([2.5 2.5],ylim,'color','k','linestyle',':');
line([8.5 8.5],ylim,'color','k','linestyle',':');
line([15.5 15.5],ylim,'color','k','linestyle',':');
title('Within sessions');
% Stats
stage_sessions = {[1,2],[3:8],[9:15],[16:21]};
y_pos = 0.35;
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
% across sessions
figure; hold on; set(gcf,'color','w','position',[200,200,300,200])
temp_var = Corr_Reward_across;
temp_var_1 = temp_var(logical(eOPN3_index_good),:);
temp_var_2 = temp_var(logical(mCherry_index_good),:);
temp_11 = nanmean(temp_var_1);
plot(temp_11,'r','linewidth',1);
temp_12 = nanstd(temp_var_1,1)./sqrt(sum(~isnan(temp_var_1)));
for ii = 1:20
    line([ii,ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','r','linewidth',1);
end
temp_21 = nanmean(temp_var_2);
plot(temp_21,'k','linewidth',1);
temp_22 = nanstd(temp_var_2,1)./sqrt(sum(~isnan(temp_var_2)));
for ii = 1:20
    line([ii,ii],[temp_21(ii)-temp_22(ii),temp_21(ii)+temp_22(ii)],'color','k','linewidth',1);
end
xlim([0.5 20.5]); xticks([1,7,13,19]); xticklabels({'1-2','7-8','13-14','19-20'}); xlabel('Session');
ylim([0 0.4]); yticks([0:0.1:0.4]); ylabel('Lever correlation');
line([1.5 1.5],ylim,'color','k','linestyle',':');
line([7.5 7.5],ylim,'color','k','linestyle',':');
line([14.5 14.5],ylim,'color','k','linestyle',':');
title('Across sessions');
% Stats
stage_sessions = {[1],[2:7],[8:14],[15:20]};
y_pos = 0.35;
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

%% Plot Figure 6J, AchR Ago injection in bilateral M1 in naive animals
load('Figure6J.mat');
% Correct rate
figure; hold on; set(gcf,'color','w','position',[200,200,200,200])
temp_var = CR(:,1:2);
temp_var_1 = temp_var(logical(Sal_index_good),:);
temp_var_2 = temp_var(logical(Ago_index_good),:);
plot(1+((rand(1,10)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,10)-0.5)*0.25),temp_var_2(:),'color',[1,0.5,0.5],'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor','k','FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor','r','FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color','k','linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color','r','linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'Sal.','Ago.'});
yticks([0:0.2:1]); yticklabels([0:20:100]); ylabel('% correct');
title('Rewarded trials');
% Stats
kk = 2;
Drug = repmat(Ago_index_good',1,kk);
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

% Lever correlation within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,200,200])
temp_var = Corr_Reward_within(:,1:2);
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(Sal_index_good),:);
temp_var_2 = temp_var(logical(Ago_index_good),:);
plot(1+((rand(1,10)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,10)-0.5)*0.25),temp_var_2(:),'color',[1,0.5,0.5],'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor','k','FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor','r','FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color','k','linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color','r','linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'Sal.','Ago.'});
ylim([-0.05 0.32]); yticks([0:0.1:0.3]); ylabel('Lever correlation');
title('Within sessions');
% Stats
kk = 2;
Drug = repmat(Ago_index_good',1,kk);
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
line([0.9,2.1],[0.3 0.3],'color','k')
text(1.5,0.31,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

% Time from cue to rewarded movement onset
figure; hold on; set(gcf,'color','w','position',[200,200,200,200])
temp_var = Cue_to_CuedRewardedMov(:,1:2);
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(Sal_index_good),:);
temp_var_2 = temp_var(logical(Ago_index_good),:);
plot(1+((rand(1,10)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,10)-0.5)*0.25),temp_var_2(:),'color',[1,0.5,0.5],'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor','k','FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor','r','FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color','k','linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color','r','linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'Sal.','Ago.'});
ylim([0 5]); yticks([0:1:5]); ylabel('Time (s)');
title('Cue to movement onset');
% Stats
kk = 2;
Drug = repmat(Ago_index_good',1,kk);
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
line([0.9,2.1],[4.75 4.75],'color','k')
text(1.5,4.85,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

% Time from rewarded movement onset to reward
figure; hold on; set(gcf,'color','w','position',[200,200,200,200])
temp_var = CuedRewardedMov_to_Reward(:,1:2);
temp_var(Trial_Num<3) = nan;
temp_var_1 = temp_var(logical(Sal_index_good),:);
temp_var_2 = temp_var(logical(Ago_index_good),:);
plot(1+((rand(1,10)-0.5)*0.25),temp_var_1(:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
plot(2+((rand(1,10)-0.5)*0.25),temp_var_2(:),'color',[1,0.5,0.5],'marker','o','linestyle','none');
temp_var_11 = nanmean(temp_var_1(:),1);
temp_var_12 = nanstd(temp_var_1(:),[],1)./sqrt(sum(~isnan(temp_var_1(:))));
temp_var_21 = nanmean(temp_var_2(:),1);
temp_var_22 = nanstd(temp_var_2(:),[],1)./sqrt(sum(~isnan(temp_var_2(:))));
bar([1],temp_var_11,'EdgeColor','k','FaceColor','none','linewidth',1,'barwidth',0.7);
bar([2],temp_var_21,'EdgeColor','r','FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_var_11-temp_var_12,temp_var_11+temp_var_12],'color','k','linewidth',1);
line([2,2],[temp_var_21-temp_var_22,temp_var_21+temp_var_22],'color','r','linewidth',1);
xlim([0.4 2.6]); xticks([1:2]); xticklabels({'Sal.','Ago.'});
ylim([0 1]); yticks([0:0.2:1]); ylabel('Time (s)');
title('Movement onset to reward');
% Stats
kk = 2;
Drug = repmat(Ago_index_good',1,kk);
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
