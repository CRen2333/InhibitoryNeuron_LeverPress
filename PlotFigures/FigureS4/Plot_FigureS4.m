%% VIP
clear
load('FigureS4.mat','VIP','-mat');
Animals = VIP.Animals;
All_Delta_df_against_N_CellTag = VIP.All_Delta_df_against_N_CellTag;
Fraction_CellTag_combo = VIP.Fraction_CellTag_combo;
Fraction_mvm_in_celltage_dec = VIP.Fraction_mvm_in_celltage_dec;
combo_Fraction_mvm_in_celltage_dec = VIP.combo_Fraction_mvm_in_celltage_dec;

%% Plot Figure S4A, significantly changed fraction
for ii = 1:4
    temp_bar(ii,:) = [sum(All_Delta_df_against_N_CellTag(:,ii)==1),sum(All_Delta_df_against_N_CellTag(:,ii)==-1),sum(All_Delta_df_against_N_CellTag(:,ii)==0)];
end
figure; hold on; set(gcf,'color','w','position',[400 400 200 200]);
bar(temp_bar(2:4,:)/sum(temp_bar(1,:))*100,'stacked','barwidth',0.7);
xlim([0.5,3.5]);ylim([0 100]); % [1,0.7,0.1;0.3,0.75,0.93;0.8,0.8,0.8]
xticks([1:3]);xticklabels(Stages(2:4));ylabel('Fraction (%)');
axis square;
% Stats
clear pval
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('M',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,3)];
[tbl,~,pval(1),~] = crosstab(x1,x2);
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(2),~] = crosstab(x1,x2);
x1 = [repmat('M',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,3);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(3),~] = crosstab(x1,x2);
[FDR] = mafdr(pval,'BHFDR', true);

%% Plot Figure S4B, significantly changed fraction, randomly select 50% animals every sweep, 
figure; hold on; set(gcf,'color','w','position',[50 50 800 200]);
for ii = 2:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_CellTag_combo{ii};
    plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
    plot(temp_bar(ii,1)/length(All_Delta_df_against_N_CellTag)*100,temp_bar(ii,2)/length(All_Delta_df_against_N_CellTag)*100,'marker','x','color','k');
    xlim([0 60]);ylim([0 80]);
    plot([0,60],[0,60],'linestyle',':','color','k');
    xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
    xlim([0 60]);xticks([0:20:60]);
    title(Stages{ii});
    axis square;
end

%% Plot Figure S4C, relationship between significalty decreased activity and movement modulation
figure; set(gcf,'color','w','position',[2000 100 600 200]); hold on
for ii = 1:3
    subplot(1,3,ii); hold on
    temp = combo_Fraction_mvm_in_celltage_dec{ii+1};
    plot(temp','color',[0.8,0.8,0.8]);
    temp = Fraction_mvm_in_celltage_dec(ii+1,:);
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
end 

%% SOM
clear
load('FigureS4.mat','SOM','-mat');
Animals = SOM.Animals;
All_Delta_df_against_N_CellTag = SOM.All_Delta_df_against_N_CellTag;
Fraction_CellTag_combo = SOM.Fraction_CellTag_combo;
Fraction_mvm_in_celltage_inc = SOM.Fraction_mvm_in_celltage_inc;
combo_Fraction_mvm_in_celltage_inc = SOM.combo_Fraction_mvm_in_celltage_inc;

%% Plot Figure S4D, significantly changed fraction
for ii = 1:4
    temp_bar(ii,:) = [sum(All_Delta_df_against_N_CellTag(:,ii)==1),sum(All_Delta_df_against_N_CellTag(:,ii)==-1),sum(All_Delta_df_against_N_CellTag(:,ii)==0)];
end
figure; hold on; set(gcf,'color','w','position',[400 400 200 200]);
bar(temp_bar(2:4,:)/sum(temp_bar(1,:))*100,'stacked','barwidth',0.7);
xlim([0.5,3.5]);ylim([0 100]); % [1,0.7,0.1;0.3,0.75,0.93;0.8,0.8,0.8]
xticks([1:3]);xticklabels(Stages(2:4));ylabel('Fraction (%)');
axis square;
% Stats
clear pval
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('M',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,3)];
[tbl,~,pval(1),~] = crosstab(x1,x2);
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(2),~] = crosstab(x1,x2);
x1 = [repmat('M',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,3);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(3),~] = crosstab(x1,x2);
[FDR] = mafdr(pval,'BHFDR', true);

%% Plot Figure S4E, significantly changed fraction, randomly select 50% animals every sweep, 
figure; hold on; set(gcf,'color','w','position',[50 50 800 200]);
for ii = 2:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_CellTag_combo{ii};
    plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
    plot(temp_bar(ii,1)/length(All_Delta_df_against_N_CellTag)*100,temp_bar(ii,2)/length(All_Delta_df_against_N_CellTag)*100,'marker','x','color','k');
    xlim([0 50]);ylim([0 50]);
    plot([0,50],[0,50],'linestyle',':','color','k');
    xlabel('Inc. fraction (%)');ylabel('Dec. fraction (%)'); axis square;
    xticks([0:25:50]);
    title(Stages{ii});
    axis square;
end

%%  Plot Figure S4F, relationship between significalty increased activity and movement modulation
figure; set(gcf,'color','w','position',[2000 100 600 200]); hold on
for ii = 1:3
    subplot(1,3,ii); hold on
    temp = combo_Fraction_mvm_in_celltage_inc{ii+1};
    plot(temp','color',[0.8,0.8,0.8]);
    temp = Fraction_mvm_in_celltage_inc(ii+1,:);
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
end 


