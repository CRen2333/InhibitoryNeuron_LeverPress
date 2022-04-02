%% Plot Figure S9A and S9B, VIP, naive stage + AChR antagonists (n&m)
load('FigureS9AS9B.mat');
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
plot(sum(All_Delta_df_against_C_CellTag(:,2)==1)/length(All_Delta_df_against_C_CellTag)*100,sum(All_Delta_df_against_C_CellTag(:,2)==-1)/length(All_Delta_df_against_C_CellTag)*100,'marker','x','color','k');
xlim([0 30]);ylim([0 30]);
plot([0,30],[0,30],'linestyle',':','color','k');
xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
xticks([0:10:30]);
yticks([0:10:30]);    
axis square;

figure; set(gcf,'color','w','position',[2000 100 200 200]); hold on
ii = 2;
temp = combo_Fraction_mvm_in_celltage_dec{ii};
plot(temp','color',[0.8,0.8,0.8]);
temp = Fraction_mvm_in_celltage_dec(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
yticks([0:0.2:1]);yticklabels([0:20:100]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction (%)');

%% Plot Figure S9C and S9D, SOM, naive stage + AChR antagonists (n&m)
load('FigureS9CS9D.mat');
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
plot(sum(All_Delta_df_against_C_CellTag(:,2)==1)/length(All_Delta_df_against_C_CellTag)*100,sum(All_Delta_df_against_C_CellTag(:,2)==-1)/length(All_Delta_df_against_C_CellTag)*100,'marker','x','color','k');
xlim([0 30]);ylim([0 30]);
plot([0,30],[0,30],'linestyle',':','color','k');
xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
xticks([0:10:30]);
yticks([0:10:30]);    
axis square;

figure; set(gcf,'color','w','position',[2000 100 200 200]); hold on
ii = 2;
temp = combo_Fraction_mvm_in_celltage_inc{ii};
plot(temp','color',[0.8,0.8,0.8]);
temp = Fraction_mvm_in_celltage_inc(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
yticks([0:0.2:1]);yticklabels([0:20:100]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction (%)');

%% Plot Figure S9E and S9F, VIP, expert stage + AChR agonists (n&m)
load('FigureS9ES9F.mat');
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
plot(sum(All_Delta_df_against_C_CellTag(:,2)==1)/length(All_Delta_df_against_C_CellTag)*100,sum(All_Delta_df_against_C_CellTag(:,2)==-1)/length(All_Delta_df_against_C_CellTag)*100,'marker','x','color','k');
xlim([0 20]);ylim([0 20]);
plot([0,20],[0,20],'linestyle',':','color','k');
xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
xticks([0:10:20]);
yticks([0:10:20]);    
axis square;

figure; set(gcf,'color','w','position',[2000 100 200 200]); hold on
ii = 2;
temp = combo_Fraction_mvm_in_celltage_inc{ii};
plot(temp','color',[0.8,0.8,0.8]);
temp = Fraction_mvm_in_celltage_inc(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
yticks([0:0.2:1]);yticklabels([0:20:100]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction (%)');

%% Plot Figure S9G and S9H, SOM, expert stage + AChR agonists (n&m)
load('FigureS9GS9H.mat');
figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
plot(sum(All_Delta_df_against_C_CellTag(:,2)==1)/length(All_Delta_df_against_C_CellTag)*100,sum(All_Delta_df_against_C_CellTag(:,2)==-1)/length(All_Delta_df_against_C_CellTag)*100,'marker','x','color','k');
xlim([0 40]);ylim([0 40]);
plot([0,40],[0,40],'linestyle',':','color','k');
xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
xticks([0:10:40]);
yticks([0:10:40]);    
axis square;

figure; set(gcf,'color','w','position',[2000 100 200 200]); hold on
ii = 2;
temp = combo_Fraction_mvm_in_celltage_dec{ii};
plot(temp','color',[0.8,0.8,0.8]);
temp = Fraction_mvm_in_celltage_dec(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
yticks([0:0.2:1]);yticklabels([0:20:100]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction (%)');

%% Plot Figure S9I, behavior at naive stage, AChR antagonists
% Correct rate
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.CR;SOM.CR];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 1]); yticks([0:0.2:1]); yticklabels([0:20:100]);
ylabel('% correct'); title('Rewarded trials');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Lever correlation, within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.Corr_within;SOM.Corr_within];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
ylabel('Lever correlation'); title('Winthin sessions');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1); % Drug or not
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,0.4,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Time from cue to movement onset
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.C_CRM;SOM.C_CRM];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 5]); yticks([0:1:5]);
ylabel('Time (s)'); title('Cue to movement onset');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,5,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Time movement onset to reward
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.RwdMVM_Rwd;SOM.RwdMVM_Rwd];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 1.2]); yticks([0:0.4:1.2]);
ylabel('Time (s)'); title('Movement onset to reward');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1.2,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

%% Plot Figure S9J, behavior at expert stage, AChR agonists
% Correct rate
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.CR;SOM.CR];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 1]); yticks([0:0.2:1]); yticklabels([0:20:100]);
ylabel('% correct'); title('Rewarded trials');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Lever correlation, within sessions
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.Corr_within;SOM.Corr_within];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
ylabel('Lever correlation'); title('Winthin sessions');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1); % Drug or not
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Time from cue to movement onset
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.C_CRM;SOM.C_CRM];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 1.2]); yticks([0:0.5:1]);
ylabel('Time (s)'); title('Cue to movement onset');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1.2,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

% Time movement onset to reward
figure; hold on; set(gcf,'color','w','position',[200,200,210,200])
temp_var = [VIP.RwdMVM_Rwd;SOM.RwdMVM_Rwd];
for ii = 1:size(temp_var,1)
    plot(temp_var(ii,:),'color',[0.8,0.8 0.8]);
end
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:2
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 2.5]);xticks([1,2]);xticklabels({'CTRL','Ago.'});
set(gca,'XColor','k','YColor','k');
ylim([0 1]); yticks([0:0.2:1]);
ylabel('Time (s)'); title('Movement onset to reward');
set(gca,'XColor','k','YColor','k'); axis square;
% stats
Sessions = repmat([1:2],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,2);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Sessions = nominal(tbl.Sessions);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(1.5,1,['p = ' num2str(pValue(2),'%.4f')],'HorizontalAlignment','center');

