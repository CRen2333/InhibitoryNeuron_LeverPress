%% Plot Figure S3, peak and duration of activity during movements
load('FigureS3.mat');
ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
for ii = 1:length(ROI)
    temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
    reordered_module(ii) = find(temp==1);
    clear temp
end
position.Naive = [1:5:5*length(ROI)];
position.Early = [2:5:5*length(ROI)];
position.Middle = [3:5:5*length(ROI)];
position.Late = [4:5:5*length(ROI)];    

%% VIP
IN = 'VIP';
PeakValue_Stage_mean = VIP.PeakValue_Stage_mean;
Duration_Stage_median = VIP.Duration_Stage_median;
color_value = VIP.color_value;
% Peak
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = PeakValue_Stage_mean{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,0.03])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Peak df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity Peak'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = PeakValue_Stage_mean{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),0.028,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.028,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.028,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.028,'*','HorizontalAlignment','center');
    end
end
% Duration
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = Duration_Stage_median{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,1.2])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('df/f FWHM','FontSize',10)
box off
title([IN ' post-movement-onset activity Duration'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = Duration_Stage_median{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),1.18,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),1.18,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),1.18,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),1.18,'*','HorizontalAlignment','center');
    end
end

%% SOM
IN = 'SOM';
PeakValue_Stage_mean = SOM.PeakValue_Stage_mean;
Duration_Stage_median = SOM.Duration_Stage_median;
color_value = SOM.color_value;
% Peak
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = PeakValue_Stage_mean{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,0.05])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Peak df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity Peak'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = PeakValue_Stage_mean{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),0.045,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.045,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.045,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.045,'*','HorizontalAlignment','center');
    end
end
% Duration
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = Duration_Stage_median{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,1.2])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('df/f FWHM','FontSize',10)
box off
title([IN ' post-movement-onset activity Duration'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = Duration_Stage_median{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),1.18,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),1.18,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),1.18,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),1.18,'*','HorizontalAlignment','center');
    end
end

%% PV
IN = 'PV';
PeakValue_Stage_mean = PV.PeakValue_Stage_mean;
Duration_Stage_median = PV.Duration_Stage_median;
color_value = PV.color_value;
% Peak
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = PeakValue_Stage_mean{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,0.045])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Peak df/f','FontSize',10)
box off
title([IN ' post-movement-onset activity Peak'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = PeakValue_Stage_mean{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),0.042,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),0.042,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),0.042,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),0.042,'*','HorizontalAlignment','center');
    end
end
% Duration
clear temp_matrix temp_mean temp_std
for ii = 1:4
    for roi = 1:16
        temp_matrix{ii}(:,roi) = Duration_Stage_median{roi}(:,ii);
    end
    temp_mean(ii,:) = nanmean(temp_matrix{ii});
    temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
end        
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
for curr_ROI = 1:length(ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
end
xlim([-1,position.Late(end)+2]); ylim([0,1.2])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('df/f FWHM','FontSize',10)
box off
title([IN ' post-movement-onset activity Duration'],'FontSize',10)    
% Statistics, REM
for curr_ROI = 1:length(ROI)
    temp_var = Duration_Stage_median{curr_ROI};
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
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
        text(position.P_ROI(curr_ROI),1.18,'****','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.001
        text(position.P_ROI(curr_ROI),1.18,'***','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.01
        text(position.P_ROI(curr_ROI),1.18,'**','HorizontalAlignment','center');
    elseif FDR(:,curr_ROI) < 0.05
        text(position.P_ROI(curr_ROI),1.18,'*','HorizontalAlignment','center');
    end
end
    