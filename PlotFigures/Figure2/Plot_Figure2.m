%% Plot Figure 2B, average traces of r-M1 and r-S1HL
load('Figure2B2C.mat');
% VIP
IN = 'VIP';
Animals = VIP.Animals;
Cued_subtract_MovActivity_all = VIP.Cued_subtract_MovActivity_all;
color_value = VIP.color_value;

figure
set(gcf,'position',[50,50,350,200]); hold on;
fields = {'Naive','Early','Middle','Late'};
subplot(1,2,1); hold on;
curr_ROI = 5; % right M1
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.002, 0.017]); ylabel('VIP mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01],'YTickLabel',{'0','0.01'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

subplot(1,2,2); hold on;
curr_ROI = 9; % righ S1HL
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.002, 0.017]); ylabel('VIP mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01],'YTickLabel',{'0','0.01'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

% SOM
IN = 'SOM';
Animals = SOM.Animals;
Cued_subtract_MovActivity_all = SOM.Cued_subtract_MovActivity_all;
color_value = SOM.color_value;

figure
set(gcf,'position',[50,50,350,200]); hold on;
fields = {'Naive','Early','Middle','Late'};
subplot(1,2,1); hold on;
curr_ROI = 5; % right M1
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.007, 0.029]); ylabel('SOM mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01,0.02],'YTickLabel',{'0','0.01','0.02'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

subplot(1,2,2); hold on;
curr_ROI = 9; % right S1HL
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.007, 0.029]); ylabel('SOM mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01,0.02],'YTickLabel',{'0','0.01','0.02'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

% PV
IN = 'PV';
Animals = PV.Animals;
Cued_subtract_MovActivity_all = PV.Cued_subtract_MovActivity_all;
color_value = PV.color_value;

figure
set(gcf,'position',[50,50,350,200]); hold on;
fields = {'Naive','Early','Middle','Late'};

subplot(1,2,1); hold on;
curr_ROI = 5; % right M1
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.004, 0.026]); ylabel('PV mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01,0.02],'YTickLabel',{'0','0.01','0.02'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

subplot(1,2,2); hold on;
curr_ROI = 9;
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',color_value(curr_field,:),'LineWidth',1);
end
xlim([1 76]); ylim([-0.004, 0.026]); ylabel('PV mean df/f'); xlabel('Time (s)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01,0.02],'YTickLabel',{'0','0.01','0.02'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

%% Plot Figure 2C, activity level in each cortical region
load('Figure2B2C.mat');
position.Naive = [1:5:5*length(ROI)];
position.Early = [2:5:5*length(ROI)];
position.Middle = [3:5:5*length(ROI)];
position.Late = [4:5:5*length(ROI)];
% VIP
IN = 'VIP';
Animals = VIP.Animals;
Cued_subtract_postActivity_all = VIP.Cued_subtract_postActivity_all;
color_value = VIP.color_value;
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module)),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module)),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module)),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module)),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5)))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5)))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5)))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5)))],'Color',color_value(4,:),'LineWidth',1);
end
xlim([-1,position.Late(end)+2]); ylim([-0.001,0.014])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Mean df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity'],'FontSize',10)

% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = [Cued_subtract_postActivity_all.Naive_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Early_AVE(:,curr_ROI),...
        Cued_subtract_postActivity_all.Middle_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Late_AVE(:,curr_ROI)];
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    p_rem_ROI(:,curr_ROI) = pValue(2);
    beta_ROI(:,curr_ROI) = beta(2);
end
beta_ROI_reorder = beta_ROI(reordered_module);
p_rem_ROI_reorder = p_rem_ROI(reordered_module);
[FDR] = mafdr(p_rem_ROI_reorder,'BHFDR', true);
position.P_ROI = (position.Early + position.Middle)/2;
hold on
for curr_ROI = 1:length(ROI)
    if FDR(:,curr_ROI) < 0.0001 
        text(position.P_ROI(curr_ROI),0.013,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.013,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.013,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.013,'*','HorizontalAlignment','center');
    end
end

% SOM
IN = 'SOM';
Animals = SOM.Animals;
Cued_subtract_postActivity_all = SOM.Cued_subtract_postActivity_all;
color_value = SOM.color_value;
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module)),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module)),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module)),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module)),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5)))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5)))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5)))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5)))],'Color',color_value(4,:),'LineWidth',1);
end
xlim([-1,position.Late(end)+2]); ylim([-0.001,0.02])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Mean df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity'],'FontSize',10)

% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = [Cued_subtract_postActivity_all.Naive_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Early_AVE(:,curr_ROI),...
        Cued_subtract_postActivity_all.Middle_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Late_AVE(:,curr_ROI)];
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    p_rem_ROI(:,curr_ROI) = pValue(2);
end
p_rem_ROI_reorder = p_rem_ROI(reordered_module);
[FDR] = mafdr(p_rem_ROI_reorder,'BHFDR', true);
position.P_ROI = (position.Early + position.Middle)/2;
hold on
for curr_ROI = 1:length(ROI)
    if FDR(:,curr_ROI) < 0.0001 
        text(position.P_ROI(curr_ROI),0.018,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.018,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.018,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.018,'*','HorizontalAlignment','center');
    end
end

% PV
IN = 'PV';
Animals = PV.Animals;
Cued_subtract_postActivity_all = PV.Cued_subtract_postActivity_all;
color_value = PV.color_value;
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module)),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module)),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module)),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module)),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Naive_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Naive_AVE(:,1)))^0.5)))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Early_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Early_AVE(:,1)))^0.5)))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Middle_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Middle_AVE(:,1)))^0.5)))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[(nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))-(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5))),...
        (nanmean(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))+(nanstd(Cued_subtract_postActivity_all.Late_AVE(:,reordered_module(curr_ROI)))./(sum(~isnan(Cued_subtract_postActivity_all.Late_AVE(:,1)))^0.5)))],'Color',color_value(4,:),'LineWidth',1);
end
xlim([-1,position.Late(end)+2]); ylim([-0.001,0.018])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Mean df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity'],'FontSize',10)

% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = [Cued_subtract_postActivity_all.Naive_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Early_AVE(:,curr_ROI),...
        Cued_subtract_postActivity_all.Middle_AVE(:,curr_ROI),Cued_subtract_postActivity_all.Late_AVE(:,curr_ROI)];
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    p_rem_ROI(:,curr_ROI) = pValue(2);
end
p_rem_ROI_reorder = p_rem_ROI(reordered_module);
[FDR] = mafdr(p_rem_ROI_reorder,'BHFDR', true);

position.P_ROI = (position.Early + position.Middle)/2;
hold on
for curr_ROI = 1:length(ROI)
    if FDR(:,curr_ROI) < 0.0001 
        text(position.P_ROI(curr_ROI),0.015,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.015,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.015,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.015,'*','HorizontalAlignment','center');
    end
end

