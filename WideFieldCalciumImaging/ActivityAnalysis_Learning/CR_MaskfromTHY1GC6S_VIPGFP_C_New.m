%% VIP_GFP control, Compare to VIP_GC6f, after ICA, no BV
clear
close all
clc

ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
for ii = 1:length(ROI)
    temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
    reordered_module(ii) = find(temp==1);
    clear temp
end

% VIP_GC6f
load('Z:\People\Chi\WFLP_IN\VIP\Craniotomy\GAP500LEEWAY150_THY1MASK\VIP_Activity_THY1MASK_OPExclude.mat', 'Cued_subtract_postActivity_all')
load('Z:\People\Chi\WFLP_IN\VIP\Craniotomy\GAP500LEEWAY150_THY1MASK\VIP_Activity_THY1MASK_OPExclude.mat', 'Cued_subtract_prepActivity_all')

fields = {'Naive','Early','Middle','Late'};
for curr_ROI = 1:length(ROI)
    for curr_field = 1:length(fields)
        field = fields{curr_field};
        Cued_subtract_MovActivity_all.(field){curr_ROI} = [Cued_subtract_prepActivity_all.(field){curr_ROI},Cued_subtract_postActivity_all.(field){curr_ROI}];
    end
end

fields = {'Naive_AVE','Early_AVE','Middle_AVE','Late_AVE'};
for ii = 1:4
    field = fields{ii};
    VIP.PostActivity{ii} = Cued_subtract_postActivity_all.(field);
end

fields = {'Naive','Early','Middle','Late'};
for ii = 1:4
    field = fields{ii};
    for curr_ROI = 1:16
        VIP.Cued_subtract_ROI_df_f_averageTrace{curr_ROI}{ii} = Cued_subtract_MovActivity_all.(field){curr_ROI};
    end
end

clear Cued_subtract_MovActivity_all Cued_subtract_postActivity_all Cued_subtract_prepActivity_all

% VIP_GFP

Initial = 'CR';
Animals = {'4383182-O','4383182-L','4383183-O'};
IN = 'VIP_GFP';

baseline_frame = [5:9];
for curr_animal = 1:3
    Animal = Animals{curr_animal};
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    data = load([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude'],'df_f_ROI_trial_arranged');
    for curr_ROI = 1:length(ROI)
        for curr_session = 1:2
            temp_matrix = data.df_f_ROI_trial_arranged{curr_ROI}{curr_session};
            baseline_mean = nanmean(temp_matrix(:,baseline_frame),2);
            temp_matrix = temp_matrix - repmat(baseline_mean,1,76);
            Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nanmean(temp_matrix,1);
        end
        VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,:) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:),1);
        VIP_GFP.session_PostActivity([curr_animal*2-1:curr_animal*2],curr_ROI) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(:,16:end),2);
        VIP_GFP.PostActivity(curr_animal,curr_ROI) = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,16:end),2);
    end
    clear data
end
    
VIP.ROI = ROI;
VIP.reordered_module = reordered_module;
VIP.Ordered_ROI = Ordered_ROI;

position.Naive = [1:6:6*length(VIP.ROI)];
position.Early = [2:6:6*length(VIP.ROI)];
position.Middle = [3:6:6*length(VIP.ROI)];
position.Late = [4:6:6*length(VIP.ROI)];
position.GFP = [5:6:6*length(VIP.ROI)];
[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];

figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,nanmean(VIP.PostActivity{1}(:,reordered_module)),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,nanmean(VIP.PostActivity{2}(:,reordered_module)),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,nanmean(VIP.PostActivity{3}(:,reordered_module)),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,nanmean(VIP.PostActivity{4}(:,reordered_module)),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
bar(position.GFP,nanmean(VIP_GFP.PostActivity(:,reordered_module)),0.20,'FaceColor',[0.5,0.5,0.5],'LineStyle','none')
for curr_ROI = 1:length(VIP.ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[(nanmean(VIP.PostActivity{1}(:,reordered_module(curr_ROI)))-(nanstd(VIP.PostActivity{1}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{1}(:,1)))^0.5))),...
        (nanmean(VIP.PostActivity{1}(:,reordered_module(curr_ROI)))+(nanstd(VIP.PostActivity{1}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{1}(:,1)))^0.5)))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[(nanmean(VIP.PostActivity{2}(:,reordered_module(curr_ROI)))-(nanstd(VIP.PostActivity{2}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{2}(:,1)))^0.5))),...
        (nanmean(VIP.PostActivity{2}(:,reordered_module(curr_ROI)))+(nanstd(VIP.PostActivity{2}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{2}(:,1)))^0.5)))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[(nanmean(VIP.PostActivity{3}(:,reordered_module(curr_ROI)))-(nanstd(VIP.PostActivity{3}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{3}(:,1)))^0.5))),...
        (nanmean(VIP.PostActivity{3}(:,reordered_module(curr_ROI)))+(nanstd(VIP.PostActivity{3}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{3}(:,1)))^0.5)))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[(nanmean(VIP.PostActivity{4}(:,reordered_module(curr_ROI)))-(nanstd(VIP.PostActivity{4}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{4}(:,1)))^0.5))),...
        (nanmean(VIP.PostActivity{4}(:,reordered_module(curr_ROI)))+(nanstd(VIP.PostActivity{4}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.PostActivity{4}(:,1)))^0.5)))],'Color',color_value(4,:),'LineWidth',1);
end
for curr_ROI = 1:length(VIP.ROI)
    line([position.GFP(curr_ROI),position.GFP(curr_ROI)],[(nanmean(VIP_GFP.PostActivity(:,reordered_module(curr_ROI)))-(nanstd(VIP_GFP.PostActivity(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP_GFP.PostActivity(:,1)))^0.5))),...
        (nanmean(VIP_GFP.PostActivity(:,reordered_module(curr_ROI)))+(nanstd(VIP_GFP.PostActivity(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP_GFP.PostActivity(:,1)))^0.5)))],'Color',[0.5,0.5,0.5],'LineWidth',1);
end

xlim([-1,position.GFP(end)+2]); ylim([-0.0012,0.014])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',position.Middle,'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Raw mean df/f','FontSize',10)
box off
title(['VIP post-movement-onset activity'],'FontSize',10)

cd('Z:\People\Chi\WFLP_IN\VIP_GFP\')
savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA.png']);

% Normalize to VIP_GC6f naive
naive_activity = nanmean(VIP.PostActivity{1}(:,reordered_module),1);
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar([1:3:3*length(VIP.ROI)],nanmean(VIP.PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1)),0.35,'FaceColor',color_value(1,:),'LineStyle','none')
bar([2:3:3*length(VIP.ROI)],nanmean(VIP_GFP.PostActivity(:,reordered_module)./repmat(naive_activity,3,1)),0.35,'FaceColor',[0.5,0.5,0.5],'LineStyle','none')
temp_1 = nanmean(VIP.PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1));
temp_2 = nanstd(VIP.PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1))/sum(~isnan(VIP.PostActivity{1}(:,1)))^0.5;
for curr_ROI = 1:length(VIP.ROI)
    line([1+(curr_ROI-1)*3,1+(curr_ROI-1)*3],[temp_1(curr_ROI)-temp_2(curr_ROI), temp_1(curr_ROI)+temp_2(curr_ROI)],'Color',color_value(1,:),'LineWidth',1);
end
temp_1 = nanmean(VIP_GFP.PostActivity(:,reordered_module)./repmat(naive_activity,3,1));
temp_2 = nanstd(VIP_GFP.PostActivity(:,reordered_module)./repmat(naive_activity,3,1))/sum(~isnan(VIP_GFP.PostActivity(:,1)))^0.5;
for curr_ROI = 1:length(VIP.ROI)
    line([2+(curr_ROI-1)*3,2+(curr_ROI-1)*3],[temp_1(curr_ROI)-temp_2(curr_ROI), temp_1(curr_ROI)+temp_2(curr_ROI)],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
for curr_ROI = 1:length(VIP.ROI)
    text(1.5+(curr_ROI-1)*3,1.2,[num2str(temp_1(curr_ROI)*100,'%.2f') '%'],'fontsize',8,'horizontalalignment','center'); 
end
xlim([-1,2+(curr_ROI-1)*3+2]); ylim([-0.25,1.3])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',[1.5:3:3*16],'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Norm. df/f','FontSize',10)
box off
title(['VIP norm. post-movement-onset activity'],'FontSize',10)

savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm.png']);

% Traces for each region
figure
set(gcf,'position',[50,50,1200,300]);
colorvalue = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
for curr_ROI = 1:length(VIP.ROI)
    subplot(2,8,curr_ROI); hold on;
    for ii = 1:4
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        temp_2 = nanstd(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})./sqrt(sum(~isnan(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colorvalue(ii,:),'FaceAlpha',0.3);
    end
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    temp_2 = nanstd(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})./sqrt(sum(~isnan(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
    for ii = 1:4
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        plot([1:76],temp_1,'color',colorvalue(ii,:),'LineWidth',1);
    end
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',[0.5 0.5 0.5],'LineWidth',1);
    
    title(VIP.Ordered_ROI{curr_ROI});
    xlim([1 76]); ylabel('Raw df/f'); xlabel('Time (sec)');
    set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
    ylim([-0.02, 0.025]); 
    set(gca,'YTick',[-0.02:0.01:0.02],'YTickLabel',[-0.02:0.01:0.02]);
    line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

end

savefig(gcf,['VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_postICA.fig']);
saveas(gcf,['VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_postICA.png']);

gfp_color = cbrewer('div','Spectral',11);
gfp_color = gfp_color(10,:);

% Only plot naive
figure;
set(gcf,'position',[50,50,600,600]);
colorvalue = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
for curr_ROI = 1:length(VIP.ROI)
    subplot(4,4,curr_ROI); hold on;
    for ii = 1
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        temp_2 = nanstd(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})./sqrt(sum(~isnan(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colorvalue(ii,:),'FaceAlpha',0.3);
    end
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    temp_2 = nanstd(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})./sqrt(sum(~isnan(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',gfp_color,'FaceAlpha',0.3);
        
    for ii = 1
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        plot([1:76],temp_1,'color',colorvalue(ii,:),'LineWidth',1);
    end    
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',gfp_color,'LineWidth',1);
    
    title(VIP.Ordered_ROI{curr_ROI},'fontsize',8);
    xlim([1 76]); ylabel('df/f'); xlabel('Time (sec)');
    set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
    ylim([-0.015, 0.022]); 
    set(gca,'YTick',[-0.02:0.01:0.02],'YTickLabel',[-0.02:0.01:0.02]);
    line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');
    axis square;
end

savefig(gcf,['Nai_VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_postICA.fig']);
saveas(gcf,['Nai_VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_postICA.png']);
print('Nai_VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_postICA.pdf','-dpdf','-bestfit'); pause(1);

% normalize to GC6f, average across regions and animals
VIP.PostActivity_corticalMean = nanmean(VIP.PostActivity{1},2);
VIP_GFP.PostActivity_corticalMean = nanmean(VIP_GFP.PostActivity,2);
naive_activity = nanmean(VIP.PostActivity_corticalMean);
figure('position',[100,100,200,200],'Color','w')
hold on;
bar([1],nanmean(VIP.PostActivity_corticalMean/naive_activity),0.7,'FaceColor',color_value(1,:),'LineStyle','none');
bar([2],nanmean(VIP_GFP.PostActivity_corticalMean/naive_activity),0.7,'FaceColor',gfp_color,'LineStyle','none');
temp_1 = nanmean(VIP.PostActivity_corticalMean/naive_activity);
temp_2 = nanstd(VIP.PostActivity_corticalMean/naive_activity)/sum(~isnan(VIP.PostActivity_corticalMean))^0.5;
line([1,1],[temp_1-temp_2, temp_1+temp_2],'Color',color_value(1,:),'LineWidth',1);
temp_1 = nanmean(VIP_GFP.PostActivity_corticalMean/naive_activity);
temp_2 = nanstd(VIP_GFP.PostActivity_corticalMean/naive_activity)/(sum(~isnan(VIP_GFP.PostActivity_corticalMean))^0.5);
line([2,2],[temp_1-temp_2, temp_1+temp_2],'Color',gfp_color,'LineWidth',1);
xlim([0.3,2.7]); ylim([-0 1.2])
set(gca,'FontSize',8,'XTick',[1,2],'XTickLabel',{'GCaMP6f','GFP'});
ylabel('Norm. df/f')

savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean.png']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean.pdf']);

% norm to GC6f, average across modules
VIP.PostActivity_corticalMean_modules = nanmean(VIP.PostActivity{1},1);
VIP_GFP.PostActivity_corticalMean_modules = nanmean(VIP_GFP.PostActivity,1);
naive_activity = VIP.PostActivity_corticalMean_modules;
figure('position',[100,100,200,200],'Color','w')
hold on;
bar([1],nanmean(VIP.PostActivity_corticalMean_modules./naive_activity),0.7,'FaceColor',color_value(1,:),'LineStyle','none');
bar([2],nanmean(VIP_GFP.PostActivity_corticalMean_modules./naive_activity),0.7,'FaceColor',gfp_color,'LineStyle','none');
temp_1 = nanmean(VIP.PostActivity_corticalMean_modules./naive_activity);
temp_2 = nanstd(VIP.PostActivity_corticalMean_modules./naive_activity)/(sum(~isnan(VIP.PostActivity_corticalMean_modules))^0.5);
line([1,1],[temp_1-temp_2, temp_1+temp_2],'Color',color_value(1,:),'LineWidth',1);
temp_1 = nanmean(VIP_GFP.PostActivity_corticalMean_modules./naive_activity);
temp_2 = nanstd(VIP_GFP.PostActivity_corticalMean_modules./naive_activity)/(sum(~isnan(VIP_GFP.PostActivity_corticalMean_modules))^0.5);
line([2,2],[temp_1-temp_2, temp_1+temp_2],'Color',gfp_color,'LineWidth',1);
xlim([0.3,2.7]); ylim([-0 1.2])
set(gca,'FontSize',8,'XTick',[1,2],'XTickLabel',{'GCaMP6f','GFP'});
ylabel('Norm. df/f')

savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean_modules.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean_modules.png']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_postICA_norm_corticalMean_modules.pdf']);

cd(['Z:\People\Chi\WFLP_IN\' IN filesep 'Craniotomy' filesep 'GAP500LEEWAY150_THY1MASK']);
if OP
    SaveName = [IN '_Activity_THY1MASK_OPExclude'];
else
    SaveName = [IN '_Activity_THY1MASK'];
end
save(SaveName,'-v7.3');
 
%% Plot example ICA modes
% VIP_GC6f: 3438544-R; VIP_GFP: CR_4383182-L
load('Z:\People\Chi\WFLP_IN\VIP\CR_3438544-R\EventAligned_Gap500\CR_3438544-R_RecICA_80.mat', 'Mode_Selected','sortMode_Retained');
load('Z:\People\Chi\WFLP_IN\VIP\CR_3438544-R\EventAligned_Gap500\ICA\ICA_80\CR_3438544-R_ICA_AllSession.mat', 'ModeICA');
Mode_exclude = ModeICA(:,~Mode_Selected);
figure; set(gcf,'color','w');
for ii = 1:size(Mode_exclude,2)
    subaxis(7,12,ii,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = Mode_exclude(:,ii);
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end

for ii = 1:size(sortMode_Retained,2)
    subaxis(7,12,ii+48,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = sortMode_Retained{ii};
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end
    
load('Z:\People\Chi\WFLP_IN\VIP_GFP\CR_4383182-L\EventAligned_Gap500\CR_4383182-L_RecICA_80.mat', 'Mode_Selected','sortMode_Retained');
load('Z:\People\Chi\WFLP_IN\VIP_GFP\CR_4383182-L\EventAligned_Gap500\ICA\ICA_80\CR_4383182-L_ICA_AllSession.mat', 'ModeICA');
Mode_exclude = ModeICA(:,~Mode_Selected);
figure; set(gcf,'color','w');
for ii = 1:size(Mode_exclude,2)
    subaxis(7,12,ii,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = Mode_exclude(:,ii);
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end

for ii = 1:size(sortMode_Retained,2)
    subaxis(7,12,ii+12*6,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = sortMode_Retained{ii};
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end
    
        