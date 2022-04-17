% Clean
clear all
close all
clc

Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};

IN = 'VIP';
OP = true;
%% Activity Level
ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
for ii = 1:length(ROI)
    temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
    reordered_module(ii) = find(temp==1);
    clear temp
end

baseline_frame = [5:9];
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    if OP
        load([Initial '_' Animal '_ROI_df_f_OPExclude'],'df_f_ROI_sub_1st_trial_arranged','df_f_ROI_trial_arranged','-mat');
    else
        load([Initial '_' Animal '_ROI_df_f'],'df_f_ROI_sub_1st_trial_arranged','df_f_ROI_trial_arranged','-mat');
    end
    load([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Final_Mask','-mat')
    for curr_ROI = 1:length(ROI)
        Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{curr_animal} = df_f_ROI_sub_1st_trial_arranged{curr_ROI};
        Cued_ROI_df_f_arranged{curr_ROI}{curr_animal} = df_f_ROI_trial_arranged{curr_ROI};
        NaivePre = nanmean(nanmean(Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{1}(:,1:15)));
        
        for curr_session = 1:11
            Cued_subtract_NaivePre_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session} = ...
                Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session} - NaivePre;
            baseline_mean = nanmean(Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session}(:,baseline_frame),2);
            Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session} = ...
                Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session} - repmat(baseline_mean,1,76);
        end            
        ROI_Mask{curr_ROI}{curr_animal} = Final_Mask{curr_ROI};
        ROI_Area(curr_ROI,curr_animal) = sum(sum(Final_Mask{curr_ROI}));
    end
    clear df_f_ROI_sub_1st_trial_arranged Final_Mask
end   

cd(['Z:\People\Chi\WFLP_IN\' IN filesep 'Craniotomy' filesep 'GAP500LEEWAY150_THY1MASK']);
if OP
    SaveName = [IN '_Activity_THY1MASK_OPExclude'];
else
    SaveName = [IN '_Activity_THY1MASK'];
end
save(SaveName,'Animals','Initial','IN','Cued_subtract_1st_ROI_df_f_arranged','Cued_ROI_df_f_arranged','Cued_subtract_ROI_df_f_arranged','ROI_Mask','ROI_Area','-v7.3');

%% Activity level
% Exclude animals with module < 10 pixels
index = sum(ROI_Area<10,1);
animalID = [1:length(Animals)];
animalID = animalID(~index);
for curr_ROI = 1:length(ROI)
    for curr_session = 1:11
        Mean_Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{curr_session} = [];
        Mean_Cued_ROI_df_f_arranged{curr_ROI}{curr_session} = [];
        Mean_Cued_subtract_NaivePre_ROI_df_f_arranged{curr_ROI}{curr_session} = [];
        Mean_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session} = [];
        for curr_animal = 1:length(animalID)
            temp = nanmean(Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{animalID(curr_animal)}{curr_session},1);
            Mean_Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{curr_session} = [Mean_Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{curr_session};temp];
            clear temp;
            temp = nanmean(Cued_ROI_df_f_arranged{curr_ROI}{animalID(curr_animal)}{curr_session},1);
            Mean_Cued_ROI_df_f_arranged{curr_ROI}{curr_session} = [Mean_Cued_ROI_df_f_arranged{curr_ROI}{curr_session};temp];
            clear temp;
            temp = nanmean(Cued_subtract_NaivePre_ROI_df_f_arranged{curr_ROI}{animalID(curr_animal)}{curr_session},1);
            Mean_Cued_subtract_NaivePre_ROI_df_f_arranged{curr_ROI}{curr_session} = [Mean_Cued_subtract_NaivePre_ROI_df_f_arranged{curr_ROI}{curr_session};temp];
            clear temp;
            temp = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{animalID(curr_animal)}{curr_session},1);
            Mean_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session} = [Mean_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session};temp];
            clear temp;
        end
    end
end
% Plot
colorvalue = colormap(hot);
for curr_animal = 1:length(animalID)
    figure(animalID(curr_animal)); set(gcf,'color','w','position',[50 50 1000 600]);
    for curr_ROI = 1:length(ROI)
        subplot(4,4,curr_ROI)
        for curr_session = 1:11
            hold on; plot(Mean_Cued_subtract_1st_ROI_df_f_arranged{curr_ROI}{curr_session}(curr_animal,:),'color',colorvalue(curr_session*5,:));
            xlim([1 76]);
        end
        title(ROI{curr_ROI});
    end
    if OP
        savefig(gcf,[Animals{animalID(curr_animal)} '_ModeActivity_OP.fig']); close all;
    else
        savefig(gcf,[Animals{animalID(curr_animal)} '_ModeActivity.fig']); close all;
    end
end

colorvalue = colormap(hot);
for curr_animal = 1:length(animalID)
    figure(animalID(curr_animal)); set(gcf,'color','w','position',[50 50 1000 600]);
    for curr_ROI = 1:length(ROI)
        subplot(4,4,curr_ROI)
        for curr_session = 1:11
            hold on; plot(Mean_Cued_ROI_df_f_arranged{curr_ROI}{curr_session}(curr_animal,:),'color',colorvalue(curr_session*5,:));
            xlim([1 76]);
        end
        title(ROI{curr_ROI});
    end
    if OP
        savefig(gcf,[Animals{animalID(curr_animal)} '_ModeActivity_nosub_OP.fig']); close all;
    else
        savefig(gcf,[Animals{animalID(curr_animal)} '_ModeActivity_nosub.fig']); close all;
    end
end

%% Activity level during movement, baseline subtracted
for curr_ROI = 1:length(ROI)
    for curr_session = 1:11
        Cued_subtract_prepActivity{curr_ROI}{curr_session} = [];
        for curr_animal = 1:length(animalID)
            if size(Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session},1) <= 1
                temp = nan(1,15);
            else
                temp = Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session}(:,1:15);
                temp = nanmean(temp,1);
            end
            Cued_subtract_prepActivity{curr_ROI}{curr_session} = [Cued_subtract_prepActivity{curr_ROI}{curr_session};temp];
            clear temp;
        end
    end
end

% Pool the data based on naive-, early-, middle-, and late-learning
for curr_ROI = 1:length(ROI)
    for curr_session = 1:11
        prepActivity_all{curr_ROI}(:,:,curr_session) = Cued_subtract_prepActivity{curr_ROI}{curr_session};
    end
    Cued_subtract_prepActivity_all.Naive{curr_ROI} = prepActivity_all{curr_ROI}(:,:,1);
    Cued_subtract_prepActivity_all.Early{curr_ROI} = nanmean(prepActivity_all{curr_ROI}(:,:,2:4),3);
    Cued_subtract_prepActivity_all.Middle{curr_ROI} = nanmean(prepActivity_all{curr_ROI}(:,:,5:8),3);
    Cued_subtract_prepActivity_all.Late{curr_ROI} = nanmean(prepActivity_all{curr_ROI}(:,:,9:11),3);
end
clear prepActivity_all

% Arrange post-movement onset activity
for curr_ROI = 1:length(ROI)
    for curr_session = 1:11
        Cued_subtract_postActivity{curr_ROI}{curr_session} = [];
        for curr_animal = animalID
            TrialNum(curr_animal,curr_session) = size(Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session},1);
            if size(Cued_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session},1) <= 1
                temp = nan(1,61);
            else
                temp = Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}{curr_session}(:,16:end);
                temp = nanmean(temp,1);
            end
            Cued_subtract_postActivity{curr_ROI}{curr_session} = [Cued_subtract_postActivity{curr_ROI}{curr_session};temp];
            clear temp;
        end
    end
end

% Pool the data based on naive-, early-, middle-, and late-learning
for curr_ROI = 1:length(ROI)
    for curr_session = 1:11
        postActivity_all{curr_ROI}(:,:,curr_session) = Cued_subtract_postActivity{curr_ROI}{curr_session};
    end
    Cued_subtract_postActivity_all.Naive{curr_ROI} = postActivity_all{curr_ROI}(:,:,1);
    Cued_subtract_postActivity_all.Early{curr_ROI} = nanmean(postActivity_all{curr_ROI}(:,:,2:4),3);
    Cued_subtract_postActivity_all.Middle{curr_ROI} = nanmean(postActivity_all{curr_ROI}(:,:,5:8),3);
    Cued_subtract_postActivity_all.Late{curr_ROI} = nanmean(postActivity_all{curr_ROI}(:,:,9:11),3);
end
clear postActivity_all

% Get average across frame
Cued_subtract_prepActivity_all.Naive_AVE = [];
Cued_subtract_prepActivity_all.Early_AVE = [];
Cued_subtract_prepActivity_all.Middle_AVE = [];
Cued_subtract_prepActivity_all.Late_AVE = [];
for curr_ROI = 1:length(ROI)
    Cued_subtract_prepActivity_all.Naive_AVE = [Cued_subtract_prepActivity_all.Naive_AVE,mean(Cued_subtract_prepActivity_all.Naive{curr_ROI},2)];
    Cued_subtract_prepActivity_all.Early_AVE = [Cued_subtract_prepActivity_all.Early_AVE,mean(Cued_subtract_prepActivity_all.Early{curr_ROI},2)];
    Cued_subtract_prepActivity_all.Middle_AVE = [Cued_subtract_prepActivity_all.Middle_AVE,mean(Cued_subtract_prepActivity_all.Middle{curr_ROI},2)];
    Cued_subtract_prepActivity_all.Late_AVE = [Cued_subtract_prepActivity_all.Late_AVE,mean(Cued_subtract_prepActivity_all.Late{curr_ROI},2)];
end

Cued_subtract_postActivity_all.Naive_AVE = [];
Cued_subtract_postActivity_all.Early_AVE = [];
Cued_subtract_postActivity_all.Middle_AVE = [];
Cued_subtract_postActivity_all.Late_AVE = [];
for curr_ROI = 1:length(ROI)
    Cued_subtract_postActivity_all.Naive_AVE = [Cued_subtract_postActivity_all.Naive_AVE,mean(Cued_subtract_postActivity_all.Naive{curr_ROI},2)];
    Cued_subtract_postActivity_all.Early_AVE = [Cued_subtract_postActivity_all.Early_AVE,mean(Cued_subtract_postActivity_all.Early{curr_ROI},2)];
    Cued_subtract_postActivity_all.Middle_AVE = [Cued_subtract_postActivity_all.Middle_AVE,mean(Cued_subtract_postActivity_all.Middle{curr_ROI},2)];
    Cued_subtract_postActivity_all.Late_AVE = [Cued_subtract_postActivity_all.Late_AVE,mean(Cued_subtract_postActivity_all.Late{curr_ROI},2)];
end

% Plot
position.Naive = [1:5:5*length(ROI)];
position.Early = [2:5:5*length(ROI)];
position.Middle = [3:5:5*length(ROI)];
position.Late = [4:5:5*length(ROI)];
[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];

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

if OP
    savefig(gcf,[IN '_PostMovActivity_sub_baseline_OPExclude.fig']);
    saveas(gcf,[IN '_PostMovActivity_sub_baseline_OPExclude.png']);
    saveas(gcf,[IN '_PostMovActivity_sub_baseline_OPExclude.eps']);
else
    savefig(gcf,[IN '_PostMovActivity.fig']); close all;
end
close(gcf);

cd(['Z:\People\Chi\WFLP_IN\' IN filesep 'Craniotomy' filesep 'GAP500LEEWAY150_THY1MASK']);
save(SaveName,'-append');

% Get averaged trace of each ROI from 4 stages
fields = {'Naive','Early','Middle','Late'};
for curr_ROI = 1:length(ROI)
    for curr_field = 1:length(fields)
        field = fields{curr_field};
        Cued_subtract_MovActivity_all.(field){curr_ROI} = [Cued_subtract_prepActivity_all.(field){curr_ROI},Cued_subtract_postActivity_all.(field){curr_ROI}];
    end
end

figure
set(gcf,'position',[50,50,1200,300]);
fields = {'Naive','Early','Middle','Late'};
colorvalue = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];

for curr_ROI = 1:length(ROI)
    subplot(2,8,curr_ROI); hold on;
    for curr_field = 1:length(fields)
        field = fields{curr_field};
        temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
        temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colorvalue(curr_field,:),'FaceAlpha',0.3);
        title(Ordered_ROI{curr_ROI});
    end
    for curr_field = 1:length(fields)
        field = fields{curr_field};
        temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
        plot([1:76],temp_1,'color',colorvalue(curr_field,:),'LineWidth',1);
    end
    xlim([1 76]); ylim([-0.005, 0.02]); ylabel('df/f'); xlabel('Time (sec)');
    set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
    set(gca,'YTick',[0,0.01,0.02,0.03],'YTickLabel',{'0','0.01','0.02','0.03'});
    line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');
end
Figure_path = ['Z:\People\Chi\WFLP_IN\' IN filesep 'Craniotomy' filesep 'GAP500LEEWAY150_THY1MASK'];
if OP
    savefig(gcf,[IN '_PostMovActivityTemporal_subtract_baseline_OPExclude.fig']);
    saveas(gcf,[IN '_PostMovActivityTemporal_subtract_baseline_OPExclude.png']);
    print([Figure_path filesep IN '_PostMovActivityTemporal_subtract_baseline_OPExclude.pdf'],'-dpdf','-bestfit'); pause(1);
else
    savefig(gcf,[IN '_PostMovActivityTemporal.fig']);
end
close(gcf);

% Plot right-M1 & right-S1HL
figure
set(gcf,'position',[50,50,350,200]); hold on;
fields = {'Naive','Early','Middle','Late'};
colorvalue = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];

subplot(1,2,1); hold on;
curr_ROI = 5; % right M1
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',colorvalue(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
    plot([1:76],temp_1,'color',colorvalue(curr_field,:),'LineWidth',1);
end
% stats
for kk = 1:76
    for ii = 2:4
        field = fields{ii};
        temp_diff = Cued_subtract_MovActivity_all.(field){curr_ROI}(animalID,kk)-Cued_subtract_MovActivity_all.Naive{curr_ROI}(animalID,kk);
        [bootstat,~] = bootstrp(20000,@mean,temp_diff);
        p_kk(ii-1,kk) = min(sum(bootstat<0),sum(bootstat>0))/20000;
    end
end
FDR = mafdr(p_kk(:),'BHFDR',true);
FDR = reshape(FDR,[3,76]);
for kk = 1:76
    if FDR(1,kk)<0.025
        plot(kk,0.0155,'.','color',[0.7,0.7,0.7])
    end
    if FDR(2,kk)<0.025
        plot(kk,0.016,'.','color',[0.4,0.4,0.4])
    end
    if FDR(3,kk)<0.025
        plot(kk,0.0165,'.','color',[0,0,0])
    end
end
xlim([1 76]); ylim([-0.002, 0.017]); ylabel('df/f'); xlabel('Time from mvm. onset (sec)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01],'YTickLabel',{'0','0.01'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

subplot(1,2,2); hold on;
curr_ROI = 9; % righ S1HL
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
    temp_2 = nanstd(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))./sqrt(sum(~isnan(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',colorvalue(curr_field,:),'FaceAlpha',0.3);
    title(Ordered_ROI{curr_ROI});
end
for curr_field = 1:length(fields)
    field = fields{curr_field};
    temp_1 = nanmean(Cued_subtract_MovActivity_all.(field){reordered_module(curr_ROI)}(animalID,:));
    plot([1:76],temp_1,'color',colorvalue(curr_field,:),'LineWidth',1);
end
% stats
for kk = 1:76
    for ii = 2:4
        field = fields{ii};
        temp_diff = Cued_subtract_MovActivity_all.(field){curr_ROI}(animalID,kk)-Cued_subtract_MovActivity_all.Naive{curr_ROI}(animalID,kk);
        [bootstat,~] = bootstrp(20000,@mean,temp_diff);
        p_kk(ii-1,kk) = min(sum(bootstat<0),sum(bootstat>0))/20000;
    end
end
FDR = mafdr(p_kk(:),'BHFDR',true);
FDR = reshape(FDR,[3,76]);
for kk = 1:76
    if FDR(1,kk)<0.025
        plot(kk,0.0155,'.','color',[0.7,0.7,0.7])
    end
    if FDR(2,kk)<0.025
        plot(kk,0.016,'.','color',[0.4,0.4,0.4])
    end
    if FDR(3,kk)<0.025
        plot(kk,0.0165,'.','color',[0,0,0])
    end
end
xlim([1 76]); ylim([-0.002, 0.017]); ylabel('df/f'); xlabel('Time from mvm. onset (sec)');
set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
set(gca,'YTick',[0,0.01],'YTickLabel',{'0','0.01'});
line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

if OP
    savefig(gcf,[IN '_rM1_S1HL_PostMovActivityTemporal_OPExclude.fig']);
    saveas(gcf,[IN '_rM1_S1HL_PostMovActivityTemporal_OPExclude.png']);
    saveas(gcf,[IN '_rM1_S1HL_PostMovActivityTemporal_OPExclude.pdf']);
else
    savefig(gcf,[IN '_rM1_PostMovActivityTemporal.fig']);
end
close(gcf);

cd(['Z:\People\Chi\WFLP_IN\' IN filesep 'Craniotomy' filesep 'GAP500LEEWAY150_THY1MASK']);
save(SaveName,'-append');
