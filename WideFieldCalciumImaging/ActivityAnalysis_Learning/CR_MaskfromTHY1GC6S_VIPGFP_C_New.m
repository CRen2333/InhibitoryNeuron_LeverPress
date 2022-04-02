clear all
close all
clc

Initial = 'CR';

Animals = {'4383182-O','4383182-L','4383183-O'};
IN = 'VIP_GFP';

OP = true;

% Get df/f
for curr_animal = 1:length(Animals)
    clearvars -except IN Initial Animals OP curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    % load mask
    load([General_Path filesep Initial '_' Animal '_Mask_from_Thy1GC6s.mat'],'Final_Mask_OP','Final_Mask_withBV_OP','-mat');
    
    All_file = dir([General_Path '\MovOnsetAligned']);
    All_filename = {All_file.name};
    All_filename = All_filename(3:end)';
    for curr_session = 1:length(All_filename)
        disp(num2str(curr_session));
        load([General_Path filesep '\MovOnsetAligned\' All_filename{curr_session}]);
        if isempty(AlignedIm_MovOnset.Index_Cued)
            Index_Cued_trial{curr_session} = [];
            Index_Cued_frame{curr_session} = [];
        else
            temp = repmat(AlignedIm_MovOnset.Index_Cued,1,76);
            temp = temp';
            Index_Cued_trial{curr_session} = AlignedIm_MovOnset.Index_Cued;
            Index_Cued_frame{curr_session} = temp(:);
        end
        if isempty(AlignedIm_MovOnset.Index_Catch)
            Index_Catch_trial{curr_session} = [];
            Index_Catch_frame{curr_session} = [];
        else
            temp = repmat(AlignedIm_MovOnset.Index_Catch,1,76);
            temp = temp';
            Index_Catch_trial{curr_session} = AlignedIm_MovOnset.Index_Catch;
            Index_Catch_frame{curr_session} = temp(:);
        end
        
        Final_Mask = Final_Mask_OP;
        % Get df/f for each ROI
        for curr_ROI = 1:length(Final_Mask)
            if isempty(AlignedIm_MovOnset.Index_Cued)
                df_f_ROI{curr_ROI}{curr_session} = [];
            else
                df_f_ROI{curr_ROI}{curr_session} = nanmean(AlignedIm_MovOnset.MovOnset_Conc(Final_Mask{curr_ROI},:),1);
            end
            trialnum{curr_session} = size(df_f_ROI{curr_ROI}{curr_session},2)/76;
            df_f_ROI_trial_arranged{curr_ROI}{curr_session} = reshape(df_f_ROI{curr_ROI}{curr_session},76,trialnum{curr_session})';
            index = logical(Index_Cued_trial{curr_session}.*~Index_Catch_trial{curr_session});
            if length(index) ~= trialnum{curr_session}
                disp(['Index:' num2str(length(index)) ', TrailNum: ' num2str(trialnum{curr_session})]);
                index = index(1:trialnum{curr_session});
            end
            Cued_ROI_df_f_arranged{curr_ROI}{curr_session} = df_f_ROI_trial_arranged{curr_ROI}{curr_session}(index,1:76);
            baseline_frame = [5:9];
            baseline_mean = nanmean(Cued_ROI_df_f_arranged{curr_ROI}{curr_session}(:,baseline_frame),2);
            Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session} = ...
                Cued_ROI_df_f_arranged{curr_ROI}{curr_session} - repmat(baseline_mean,1,76);
            
            Cued_ROI_df_f_average{curr_ROI}(curr_session,:) = nanmean(Cued_ROI_df_f_arranged{curr_ROI}{curr_session},1);
            Cued_subtract_ROI_df_f_average{curr_ROI}(curr_session,:) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
        end
        
        Final_Mask = Final_Mask_withBV_OP;
        % Get df/f for each ROI
        for curr_ROI = 1:length(Final_Mask)
            if isempty(AlignedIm_MovOnset.Index_Cued)
                withBV_df_f_ROI{curr_ROI}{curr_session} = [];
            else
                withBV_df_f_ROI{curr_ROI}{curr_session} = nanmean(AlignedIm_MovOnset.MovOnset_Conc(Final_Mask{curr_ROI},:),1);
            end
            trialnum{curr_session} = size(withBV_df_f_ROI{curr_ROI}{curr_session},2)/76;
            withBV_df_f_ROI_trial_arranged{curr_ROI}{curr_session} = reshape(withBV_df_f_ROI{curr_ROI}{curr_session},76,trialnum{curr_session})';
            index = logical(Index_Cued_trial{curr_session}.*~Index_Catch_trial{curr_session});
            if length(index) ~= trialnum{curr_session}
                disp(['Index:' num2str(length(index)) ', TrailNum: ' num2str(trialnum{curr_session})]);
                index = index(1:trialnum{curr_session});
            end
            withBV_Cued_ROI_df_f_arranged{curr_ROI}{curr_session} = withBV_df_f_ROI_trial_arranged{curr_ROI}{curr_session}(index,1:76);
            baseline_frame = [5:9];
            baseline_mean = nanmean(withBV_Cued_ROI_df_f_arranged{curr_ROI}{curr_session}(:,baseline_frame),2);
            withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session} = ...
                withBV_Cued_ROI_df_f_arranged{curr_ROI}{curr_session} - repmat(baseline_mean,1,76);
            
            withBV_Cued_ROI_df_f_average{curr_ROI}(curr_session,:) = nanmean(withBV_Cued_ROI_df_f_arranged{curr_ROI}{curr_session},1);
            withBV_Cued_subtract_ROI_df_f_average{curr_ROI}(curr_session,:) = nanmean(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
        end
    end  
    clear Final_Mask_OP Final_Mask_withBV_OP AlignedIm_MovOnset
    
    save([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude_NOICA'],'-v7.3');
    
end

%% Activity Level, VIP_GFP
clear all
close all
clc

ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
for ii = 1:length(ROI)
    temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
    reordered_module(ii) = find(temp==1);
    clear temp
end

Initial = 'CR';
Animals = {'4383182-O','4383182-L','4383183-O'};
IN = 'VIP_GFP';

% no BV
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    data = load([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude_NOICA'],'Cued_subtract_ROI_df_f_arranged');
    TrialNum(curr_animal,1) = size(data.Cued_subtract_ROI_df_f_arranged{1}{1},1);
    TrialNum(curr_animal,2) = size(data.Cued_subtract_ROI_df_f_arranged{1}{2},1);
    for curr_ROI = 1:length(ROI)
        for curr_session = 1:2
            Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nanmean(data.Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
        end
        Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,:) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:),1);
        session_PostActivity([curr_animal*2-1:curr_animal*2],curr_ROI) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(:,16:end),2);
        PostActivity(curr_animal,curr_ROI) = nanmean(Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,16:end),2);
    end
    clear data
end

% Animal
figure; hold on
for ii = 1:16
    subplot(4,4,ii); hold on; 
    plot(Cued_subtract_ROI_df_f_averageTrace{curr_ROI}'); title(Ordered_ROI(ii));
end

% Session
figure; hold on
for ii = 1:16
    subplot(4,4,ii); hold on;
    temp = cell2mat(Cued_subtract_ROI_df_f_arranged{curr_ROI}');
    plot(temp'); title(Ordered_ROI(ii));
end

% with BV
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    data = load([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude_NOICA'],'withBV_Cued_subtract_ROI_df_f_arranged');
    TrialNum(curr_animal,1) = size(data.withBV_Cued_subtract_ROI_df_f_arranged{1}{1},1);
    TrialNum(curr_animal,2) = size(data.withBV_Cued_subtract_ROI_df_f_arranged{1}{2},1);
    for curr_ROI = 1:length(ROI)
        for curr_session = 1:2
            withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nanmean(data.withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
        end
        withBV_Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,:) = nanmean(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:),1);
        withBV_session_PostActivity([curr_animal*2-1:curr_animal*2],curr_ROI) = nanmean(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(:,16:end),2);
        withBV_PostActivity(curr_animal,curr_ROI) = nanmean(withBV_Cued_subtract_ROI_df_f_averageTrace{curr_ROI}(curr_animal,16:end),2);
    end
    clear data
end

% Animal
figure; hold on
for ii = 1:16
    subplot(4,4,ii); hold on; 
    plot(withBV_Cued_subtract_ROI_df_f_averageTrace{curr_ROI}'); title(Ordered_ROI(ii));
end

% Session
figure; hold on
for ii = 1:16
    subplot(4,4,ii); hold on;
    temp = cell2mat(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}');
    plot(temp'); title(Ordered_ROI(ii));
end


cd('Z:\People\Chi\WFLP_IN\VIP_GFP');
save('VIP_GFP_Activity_THY1MASK_OPExclude_NOICA','-v7.3');

%% Activity level, VIP
clear all
close all
clc

ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
for ii = 1:length(ROI)
    temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
    reordered_module(ii) = find(temp==1);
    clear temp
end

% no BV
Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};
IN = 'VIP';
stage_sessions = {[1],[2:4],[5:8],[9:11]};
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    data = load([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude_NOICA'],'Cued_subtract_ROI_df_f_arranged');
    for curr_ROI = 1:length(ROI)
        for curr_session = 1:11
            TrialNum(curr_animal,curr_session) = size(data.Cued_subtract_ROI_df_f_arranged{1}{curr_session},1);
            if TrialNum(curr_animal,curr_session)<1
                Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nan;
            else
                Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nanmean(data.Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
            end
        end
        session_PostActivity{curr_ROI}(curr_animal,:) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(:,16:end),2);
        for ii_stage = 1:4
            Cued_subtract_ROI_df_f_averageTrace{curr_ROI}{ii_stage}(curr_animal,:) = nanmean(Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(stage_sessions{ii_stage},:),1);
            PostActivity{ii_stage}(curr_animal,curr_ROI) = nanmean(session_PostActivity{curr_ROI}(curr_animal,stage_sessions{ii_stage}),2);
        end
    end
    clear data
end

% with BV
Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};
IN = 'VIP';
stage_sessions = {[1],[2:4],[5:8],[9:11]};
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    General_Path = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500'];
    data = load([General_Path filesep Initial '_' Animal '_ROI_df_f_OPExclude_NOICA'],'withBV_Cued_subtract_ROI_df_f_arranged');
    for curr_ROI = 1:length(ROI)
        for curr_session = 1:11
            TrialNum(curr_animal,curr_session) = size(data.withBV_Cued_subtract_ROI_df_f_arranged{1}{curr_session},1);
            if TrialNum(curr_animal,curr_session)<1
                withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nan;
            else
                withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(curr_session,:) = nanmean(data.withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_session},1);
            end
        end
        withBV_session_PostActivity{curr_ROI}(curr_animal,:) = nanmean(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(:,16:end),2);
        for ii_stage = 1:4
            withBV_Cued_subtract_ROI_df_f_averageTrace{curr_ROI}{ii_stage}(curr_animal,:) = nanmean(withBV_Cued_subtract_ROI_df_f_arranged{curr_ROI}{curr_animal}(stage_sessions{ii_stage},:),1);
            withBV_PostActivity{ii_stage}(curr_animal,curr_ROI) = nanmean(withBV_session_PostActivity{curr_ROI}(curr_animal,stage_sessions{ii_stage}),2);
        end
    end
    clear data
end

cd('Z:\People\Chi\WFLP_IN\VIP_GFP');
save('VIP_Activity_THY1MASK_OPExclude_NOICA','-v7.3');

%% Plot
close all
clear
clc

cd('Z:\People\Chi\WFLP_IN\VIP_GFP');
VIP = load('VIP_Activity_THY1MASK_OPExclude_NOICA');
VIP_GFP = load('VIP_GFP_Activity_THY1MASK_OPExclude_NOICA');

% no BV, post activity
position.Naive = [1:6:6*length(VIP.ROI)];
position.Early = [2:6:6*length(VIP.ROI)];
position.Middle = [3:6:6*length(VIP.ROI)];
position.Late = [4:6:6*length(VIP.ROI)];
position.GFP = [5:6:6*length(VIP.ROI)];
[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
% color_value = [196,223,163;166,198,128;119,155,76;98,130,59]./255;

reordered_module = VIP.reordered_module;

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

xlim([-1,position.GFP(end)+2]); ylim([-0.001,0.015])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',position.Middle,'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Raw mean df/f','FontSize',10)
box off
title(['VIP post-movement-onset activity'],'FontSize',10)

savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_NOICA.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_NOICA.png']);

% Normalize to naive
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
xlim([-1,2+(curr_ROI-1)*3+2]); ylim([-0.15,1.3])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',[1.5:3:3*16],'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Norm. df/f','FontSize',10)
box off
title(['VIP norm. post-movement-onset activity'],'FontSize',10)

savefig(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_NOICA_norm.fig']);
saveas(gcf,['VIP_PostMovActivity_sub_baseline_OPExclude_NOICA_norm.png']);

% Traces
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

savefig(gcf,['VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_NOICA.fig']);
saveas(gcf,['VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_NOICA.png']);

% with BV, post activity
position.Naive = [1:6:6*length(VIP.ROI)];
position.Early = [2:6:6*length(VIP.ROI)];
position.Middle = [3:6:6*length(VIP.ROI)];
position.Late = [4:6:6*length(VIP.ROI)];
position.GFP = [5:6:6*length(VIP.ROI)];
[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
% color_value = [196,223,163;166,198,128;119,155,76;98,130,59]./255;

reordered_module = VIP.reordered_module;

figure('position',[100,100,1200,200],'Color','w')
hold on;
bar(position.Naive,nanmean(VIP.withBV_PostActivity{1}(:,reordered_module)),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
bar(position.Early,nanmean(VIP.withBV_PostActivity{2}(:,reordered_module)),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
bar(position.Middle,nanmean(VIP.withBV_PostActivity{3}(:,reordered_module)),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
bar(position.Late,nanmean(VIP.withBV_PostActivity{4}(:,reordered_module)),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
bar(position.GFP,nanmean(VIP_GFP.withBV_PostActivity(:,reordered_module)),0.20,'FaceColor',[0.5,0.5,0.5],'LineStyle','none')
for curr_ROI = 1:length(VIP.ROI)
    line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[(nanmean(VIP.withBV_PostActivity{1}(:,reordered_module(curr_ROI)))-(nanstd(VIP.withBV_PostActivity{1}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{1}(:,1)))^0.5))),...
        (nanmean(VIP.withBV_PostActivity{1}(:,reordered_module(curr_ROI)))+(nanstd(VIP.withBV_PostActivity{1}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{1}(:,1)))^0.5)))],'Color',color_value(1,:),'LineWidth',1);
    line([position.Early(curr_ROI),position.Early(curr_ROI)],[(nanmean(VIP.withBV_PostActivity{2}(:,reordered_module(curr_ROI)))-(nanstd(VIP.withBV_PostActivity{2}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{2}(:,1)))^0.5))),...
        (nanmean(VIP.withBV_PostActivity{2}(:,reordered_module(curr_ROI)))+(nanstd(VIP.withBV_PostActivity{2}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{2}(:,1)))^0.5)))],'Color',color_value(2,:),'LineWidth',1);
    line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[(nanmean(VIP.withBV_PostActivity{3}(:,reordered_module(curr_ROI)))-(nanstd(VIP.withBV_PostActivity{3}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{3}(:,1)))^0.5))),...
        (nanmean(VIP.withBV_PostActivity{3}(:,reordered_module(curr_ROI)))+(nanstd(VIP.withBV_PostActivity{3}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{3}(:,1)))^0.5)))],'Color',color_value(3,:),'LineWidth',1);
    line([position.Late(curr_ROI),position.Late(curr_ROI)],[(nanmean(VIP.withBV_PostActivity{4}(:,reordered_module(curr_ROI)))-(nanstd(VIP.withBV_PostActivity{4}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{4}(:,1)))^0.5))),...
        (nanmean(VIP.withBV_PostActivity{4}(:,reordered_module(curr_ROI)))+(nanstd(VIP.withBV_PostActivity{4}(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP.withBV_PostActivity{4}(:,1)))^0.5)))],'Color',color_value(4,:),'LineWidth',1);
end
for curr_ROI = 1:length(VIP.ROI)
    line([position.GFP(curr_ROI),position.GFP(curr_ROI)],[(nanmean(VIP_GFP.withBV_PostActivity(:,reordered_module(curr_ROI)))-(nanstd(VIP_GFP.withBV_PostActivity(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP_GFP.withBV_PostActivity(:,1)))^0.5))),...
        (nanmean(VIP_GFP.withBV_PostActivity(:,reordered_module(curr_ROI)))+(nanstd(VIP_GFP.withBV_PostActivity(:,reordered_module(curr_ROI)))./(sum(~isnan(VIP_GFP.withBV_PostActivity(:,1)))^0.5)))],'Color',[0.5,0.5,0.5],'LineWidth',1);
end

xlim([-1,position.GFP(end)+2]); ylim([-0.001,0.018])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',position.Middle,'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Raw mean df/f','FontSize',10)
box off
title(['VIP post-movement-onset activity'],'FontSize',10)

savefig(gcf,['withBV_VIP_PostMovActivity_sub_baseline_OPExclude_NOICA.fig']);
saveas(gcf,['withBV_VIP_PostMovActivity_sub_baseline_OPExclude_NOICA.png']);

% Normalize to naive
naive_activity = nanmean(VIP.withBV_PostActivity{1}(:,reordered_module),1);
figure('position',[100,100,1200,200],'Color','w')
hold on;
bar([1:3:3*length(VIP.ROI)],nanmean(VIP.withBV_PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1)),0.35,'FaceColor',color_value(1,:),'LineStyle','none')
bar([2:3:3*length(VIP.ROI)],nanmean(VIP_GFP.withBV_PostActivity(:,reordered_module)./repmat(naive_activity,3,1)),0.35,'FaceColor',[0.5,0.5,0.5],'LineStyle','none')
temp_1 = nanmean(VIP.withBV_PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1));
temp_2 = nanstd(VIP.withBV_PostActivity{1}(:,reordered_module)./repmat(naive_activity,11,1))/sum(~isnan(VIP.withBV_PostActivity{1}(:,1)))^0.5;
for curr_ROI = 1:length(VIP.ROI)
    line([1+(curr_ROI-1)*3,1+(curr_ROI-1)*3],[temp_1(curr_ROI)-temp_2(curr_ROI), temp_1(curr_ROI)+temp_2(curr_ROI)],'Color',color_value(1,:),'LineWidth',1);
end
temp_1 = nanmean(VIP_GFP.withBV_PostActivity(:,reordered_module)./repmat(naive_activity,3,1));
temp_2 = nanstd(VIP_GFP.withBV_PostActivity(:,reordered_module)./repmat(naive_activity,3,1))/sum(~isnan(VIP_GFP.withBV_PostActivity(:,1)))^0.5;
for curr_ROI = 1:length(VIP.ROI)
    line([2+(curr_ROI-1)*3,2+(curr_ROI-1)*3],[temp_1(curr_ROI)-temp_2(curr_ROI), temp_1(curr_ROI)+temp_2(curr_ROI)],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
for curr_ROI = 1:length(VIP.ROI)
    text(1.5+(curr_ROI-1)*3,1.2,[num2str(temp_1(curr_ROI)*100,'%.2f') '%'],'fontsize',8,'horizontalalignment','center'); 
end
xlim([-1,2+(curr_ROI-1)*3+2]); ylim([-0.15,1.3])
set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',[1.5:3:3*16],'XTickLabel',VIP.Ordered_ROI);
xlabel('Cortical module','FontSize',10);ylabel('Norm. df/f','FontSize',10)
box off
title(['VIP norm. post-movement-onset activity'],'FontSize',10)

savefig(gcf,['withBV_VIP_PostMovActivity_sub_baseline_OPExclude_NOICA_norm.fig']);
saveas(gcf,['withBV_VIP_PostMovActivity_sub_baseline_OPExclude_NOICA_norm.png']);

% Traces
figure
set(gcf,'position',[50,50,1200,300]);
colorvalue = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
for curr_ROI = 1:length(VIP.ROI)
    subplot(2,8,curr_ROI); hold on;
    for ii = 1:4
        temp_1 = nanmean(VIP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        temp_2 = nanstd(VIP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})./sqrt(sum(~isnan(VIP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colorvalue(ii,:),'FaceAlpha',0.3);
    end
    temp_1 = nanmean(VIP_GFP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    temp_2 = nanstd(VIP_GFP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})./sqrt(sum(~isnan(VIP_GFP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
    for ii = 1:4
        temp_1 = nanmean(VIP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        plot([1:76],temp_1,'color',colorvalue(ii,:),'LineWidth',1);
    end
    temp_1 = nanmean(VIP_GFP.withBV_Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',[0.5 0.5 0.5],'LineWidth',1);
    
    title(VIP.Ordered_ROI{curr_ROI});
    xlim([1 76]); ylabel('Raw df/f'); xlabel('Time (sec)');
    set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
    ylim([-0.02, 0.025]); 
    set(gca,'YTick',[-0.02:0.01:0.03],'YTickLabel',[-0.02:0.01:0.02]);
    line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');

end

savefig(gcf,['withBV_VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_NOICA.fig']);
saveas(gcf,['withBV_VIP_PostMovActivityTemporal_subtract_baseline_OPExclude_NOICA.png']);

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

% norm, average across animals
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

% norm, average across modules
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
    
        