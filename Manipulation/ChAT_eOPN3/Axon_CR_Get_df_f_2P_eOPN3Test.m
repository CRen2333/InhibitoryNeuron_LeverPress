%% Group axonal ROIs belonging to the same axon

clear clusters
Date = '210622';
cd(['Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT\CR_4333447-R\' Date '\Summed_2P_image']);
load(['CR_' Date '_4333447-R_Rec1_001_001_summed_50_template.roi'],'-mat');
clear polygon

% 0622, 10 mW, 15 sec
clusters{1} = [9:11,15]; clusters{2} = [16:20]; clusters{3} = [23,28,29]; clusters{4} = [24,25,31]; clusters{5} = [26,27];
clusters{6} = [36,37,39]; clusters{7} = [42:47]; clusters{8} = [50,52,53]; clusters{9} = [54:58];

Sessions = {'Rec_2','Rec_3','Rec_4'};

for ii_rec = 1:length(Sessions)
    load(['CR_' Date '_4333447-R_Rec' Sessions{ii_rec}(end) '_001_001_summed_50_template.roi'],'polygon','-mat');
    for ii = 1:length(clusters)
        polygon_cluster.ROI{ii} = cell2mat(polygon.ROI(clusters{ii})');
    end
    in_cluster = cell2mat(clusters);
    rest_roi = find(~ismember([1:length(polygon.ROI)],in_cluster));
    polygon_cluster.ROI = [polygon_cluster.ROI, polygon.ROI(rest_roi)];
    clear polygon
    polygon = polygon_cluster;
    save(['CR_' Date '_4333447-R_Rec' Sessions{ii_rec}(end) '_001_001_summed_50_template_cluster.mat'],'polygon','-v7.3');
    clear polygon polygon_cluster
end

%% Get roi traces & df/f
clear all;
close all;
clc;

Animal = 'CR_4333447-R';
Date = '210622';
Rig = 'MOM';
disp(Animal);
General_path = fullfile('Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT',Animal,Date);
Local_path = fullfile('F:\Data\MotionCorrection',Animal,Date);
ROI_path = ['Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT\' Animal '\' Date '\Summed_2P_image'];
disp(Date);            
Target_path = [General_path filesep 'df_f'];
if ~exist(Target_path)
    mkdir(Target_path);
end

image_folder_path = fullfile(Local_path);
tiff_folders = dir(image_folder_path);
tiff_folders = {tiff_folders.name};
tiff_folders = sort(tiff_folders);
tiff_folders = tiff_folders(3:end);
tiff_folders = {'Rec_2','Rec_3','Rec_4'};

for curr_session = 1:length(tiff_folders)
    roi_filename = fullfile(ROI_path,[Animal(1:3) Date Animal(3:end) '_Rec' tiff_folders{curr_session}(end) '_001_001_summed_50_template_cluster.mat']);
    bg_roi_filename = fullfile(ROI_path,[Animal(1:3) Date Animal(3:end) '_Rec' tiff_folders{curr_session}(end) '_001_001_summed_50_bg.roi']);
    tiff_path = fullfile(image_folder_path,tiff_folders{curr_session},'Z1','motioncorrected_tiff');
    interp_sub = 0;
    local_comp = 0;
    switch Rig
        case 'MOM'
            framerate = 28; 
            [roi_trace{curr_session}] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp);     
            [roi_trace_bg{curr_session}] = AP_getConcatTrace_continuous_batch(bg_roi_filename,tiff_path,interp_sub,local_comp);                 
        case 'BScope1'
            framerate = 30;
            [roi_trace{curr_session},roi_trace_bg{curr_session},polygon{curr_session}] = CR_Suite2P_Axon_Rig3_getConcatTrace_continuous_batch(roi_filename,bg_roi_file,tiff_path,interp_sub,local_comp);                        
    end

    % select 10 rois for qualirt check
    if size(roi_trace{curr_session},1) < 10                    
        selected_roi_index = [1:5];
    elseif size(roi_trace{curr_session},1) < 20                    
        selected_roi_index = [1:2:10];
    else
        selected_roi_index = [1:2:20];
    end
    roi_for_qualitycheck{curr_session} = roi_trace{curr_session}(selected_roi_index,:);
end
Target_filename = [Animal '_' Date '_ROI_Traces_cluster.mat'];
save([Target_path filesep Target_filename],'roi_trace','roi_trace_bg','-v7.3');
clear roi_trace roi_trace_bg 
save([General_path filesep 'df_f' filesep Animal '_' Date '_ROI_Traces_Check_cluster.mat'],'roi_for_qualitycheck','-v7.3');

%% Check truncate point 
clear all;
close all;
clc;

Animal = 'CR_4333447-R';
Date = '210622';
Rig = 'MOM';
disp(Animal);
General_path = fullfile('Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT',Animal,Date);
load([General_path filesep 'df_f' filesep Animal '_' Date '_ROI_Traces_Check_cluster.mat'],'-mat');            

for curr_session = 1:length(roi_for_qualitycheck)
    window_start = 1:1000;
    window_end = size(roi_for_qualitycheck{curr_session},2)-1000:size(roi_for_qualitycheck{curr_session},2);
    figure; set(gcf,'position',[0,0,1200,600]); hold on;
    for curr_roi = 1:size(roi_for_qualitycheck{curr_session},1)
        plot(roi_for_qualitycheck{curr_session}(curr_roi,window_start)-500*(curr_roi-1));
    end
    title([Animal ' ' num2str(curr_session) '/' num2str(length(roi_for_qualitycheck)) ' Start']);
    [temp_start,~] = ginput(1);
    temp_start = round(temp_start);
    close(gcf);
    figure; set(gcf,'position',[50,50,1400,800]); hold on;
    for curr_roi = 1:size(roi_for_qualitycheck{curr_session},1)
        plot(roi_for_qualitycheck{curr_session}(curr_roi,window_end)-500*(curr_roi-1));
    end
    title([Animal ' ' num2str(curr_session) '/' num2str(length(roi_for_qualitycheck)) ' End']);
    [temp_end,~] = ginput(1);
    temp_end = temp_end + size(roi_for_qualitycheck{curr_session},2) - 1000;
    temp_end = round(temp_end);
    close(gcf);
    figure; set(gcf,'position',[50,50,1400,800]); hold on;
    for curr_roi = 1:size(roi_for_qualitycheck{curr_session},1)
        plot(roi_for_qualitycheck{curr_session}(curr_roi,:)-1000*(curr_roi-1));
    end
    title([Animal ' ' num2str(curr_session) '/' num2str(length(roi_for_qualitycheck))]);
    [temp_all,~] = ginput(2);
    temp_all = round(temp_all);
    close(gcf)
    % start
    truncatePoint{curr_session}(1) = max([temp_start,temp_all(1),1]);
    % end
    truncatePoint{curr_session}(2) = min([temp_end,temp_all(2),size(roi_for_qualitycheck{curr_session},2)]);
    if truncatePoint{curr_session}(2) == size(roi_for_qualitycheck{curr_session},2)
        istruncate{curr_session} = true;
        truncatePoint{curr_session}(2) = truncatePoint{curr_session}(2)-20; % last 20 data points are often bad;
    else
        istruncate{curr_session} = true;
    end
end
save([General_path filesep 'df_f' filesep Animal '_' Date '_ROI_Traces_cluster.mat'],'truncatePoint','istruncate','-append');
clear roi_trace roi_trace_bg roi_for_qualitycheck truncatePoint istruncate temp_start temp_end temp_all
       
%% Estimate baseline and get df/f
clear

Animal = 'CR_4333447-R';
Date = '210622';
Rig = 'MOM';
disp(Animal);
General_path = fullfile('Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT',Animal,Date);
cd([General_path filesep 'df_f']);
load([Animal '_' Date '_ROI_Traces_cluster.mat'],'-mat');            

df_threshold = 100;
if ~exist('roi_trace_origin')
    disp('Saving original traces...');
    roi_trace_origin = roi_trace;
    roi_trace_bg_origin = roi_trace_bg;
end
for curr_session = 1:length(roi_trace) 
    disp([Animal ' Session' num2str(curr_session)]);
    roi_trace{curr_session} = roi_trace_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2));
    roi_trace_bg{curr_session} = roi_trace_bg_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2));
    roi_trace_bg{curr_session} = nanmean(roi_trace_bg{curr_session},1);
    roi_trace_long{curr_session} = roi_trace{curr_session} - repmat(roi_trace_bg{curr_session},size(roi_trace{curr_session},1),1);
    framerate = 28;
    roi_trace_df_2{curr_session} = nan(size(roi_trace_long{curr_session}));
    roi_trace_baseline_2{curr_session} = nan(size(roi_trace_long{curr_session}));
    roi_trace_df_3{curr_session} = nan(size(roi_trace_long{curr_session}));
    roi_trace_baseline_3{curr_session} = nan(size(roi_trace_long{curr_session}));
    all_nan_index = sum(isnan(roi_trace_long{curr_session}),2);
    baseline_window = round(framerate*1*60);
    prctile_thre = 30;
    tic
    for frame = 1:size(roi_trace_long{curr_session},2)
        if frame <= baseline_window
            roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,1:frame+baseline_window),prctile_thre,2);
        elseif frame >= size(roi_trace_long{curr_session},2)-baseline_window
            roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,frame-baseline_window:size(roi_trace_long{curr_session},2)),prctile_thre,2);
        else
            roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,frame-baseline_window:frame+baseline_window),prctile_thre,2);
        end
    end
    if any(all_nan_index)
        disp('Calculating for nan containing traces');
        nan_rois = find(all_nan_index);
        for roi = nan_rois
            temp_trace = roi_trace_long{curr_session}(roi,:);
            nan_index = ~isnan(temp_trace);
            if sum(~nan_index) == length(temp_trace)
                disp([num2str(roi) ' all nan, ignore']);
                continue
            end
            temp_trace = temp_trace(nan_index);
            clear temp_baseline
            for frame = 1:length(temp_trace)
                if frame <= baseline_window
                    temp_baseline(frame) = prctile(temp_trace(1:frame+baseline_window),prctile_thre);
                elseif frame >= size(roi_trace_long{curr_session},2)-baseline_window
                    temp_baseline(frame) = prctile(temp_trace(frame-baseline_window:size(roi_trace_long{curr_session},2)),prctile_thre);
                else
                    temp_baseline(frame) = prctile(temp_trace(frame-baseline_window:frame+baseline_window),prctile_thre);
                end
            end
            roi_trace_baseline_2{curr_session}(roi,nan_index) = temp_baseline;
        end
    end
    toc
    roi_trace_df_2{curr_session} = (roi_trace_long{curr_session}-roi_trace_baseline_2{curr_session})./roi_trace_baseline_2{curr_session};
    % get rid of noise
    roi_trace_df_2{curr_session}(abs(roi_trace_df_2{curr_session}) >= df_threshold) = nan;

    % Plot for post check
    FigTargetPath = [General_path filesep 'df_f' filesep 'TraceToCheck'];
    if ~exist(FigTargetPath)
        mkdir(FigTargetPath)
    end
    close all;
    figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 900]);
    hold on;
    for curr_roi = 1:size(roi_trace_long{curr_session},1)
        tt = [1:length(roi_trace_long{curr_session})]/framerate;
        plot(tt,roi_trace_df_2{curr_session}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
    end
    set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');
    savefig(gcf,[FigTargetPath filesep Animal '_' Date '_' num2str(curr_session) '_cluster.fig']);
    saveas(gcf,[FigTargetPath filesep Animal '_' Date '_' num2str(curr_session) '_cluster.png']);
    pause(0.1);
    close(gcf);

end
save([General_path filesep 'df_f' filesep Animal '_' Date '_ROI_Traces_cluster.mat'],'roi_trace_origin','roi_trace_bg_origin','roi_trace_long','roi_trace_df_2','roi_trace_baseline_2','-append'); 

%%
clear
Date = '210622';
cd(['Z:\People\Chi\TwoP_IN\VirusTest\eOPN3_ChAT\CR_4333447-R\' Date '\df_f']);
load(['CR_4333447-R_' Date '_ROI_Traces_cluster.mat'], 'roi_trace_df_2');
% Reject too dim ROIs
RejectedROI = [15,20,24];
RejectedROI_G = [12];

% Align to ca event onset based on averaged traces
clear df_f_aligned_CaOn
for curr_session = 1:length(roi_trace_df_2)
    noRejectedROI = ~ismember([1:size(roi_trace_df_2{curr_session},1)],unique([RejectedROI,RejectedROI_G]));
    temp_mean_trace = nanmean(roi_trace_df_2{curr_session}(noRejectedROI,:),1);
    roi_trace_df_2{curr_session} = roi_trace_df_2{curr_session}(noRejectedROI,:);
    threshold = 1; frame_rate = 28; smooth_window_sec = 3;
    [~,~,caEvents] = CR_GetCaEvents_ChAT_axon(temp_mean_trace',threshold,frame_rate,smooth_window_sec,false);
    caEvents = caEvents'>0;
    caEvents_on = find(diff([0,caEvents])==1);
    caEvents_on(caEvents_on<1*28)= [];
    caEvents_on(caEvents_on>120*28-3*28)= [];
    for ii = 1:length(caEvents_on)
        curr_onset = caEvents_on(ii);
        for curr_roi = 1:size(roi_trace_df_2{curr_session},1)
            df_f_aligned_CaOn{curr_session}{curr_roi}(:,ii) = roi_trace_df_2{curr_session}(curr_roi,curr_onset-28*1:curr_onset+28*3);
        end
    end
end

for curr_session = 1:length(df_f_aligned_CaOn)
    for curr_roi = 1:length(df_f_aligned_CaOn{curr_session})
        df_f_aligned_CaOn_MeanXTrial{curr_session}(curr_roi,:) = nanmean(df_f_aligned_CaOn{curr_session}{curr_roi},2);
        % subtract baseline
        df_f_aligned_CaOn_MeanXTrial_sub{curr_session}(curr_roi,:) = df_f_aligned_CaOn_MeanXTrial{curr_session}(curr_roi,:)- nanmean(df_f_aligned_CaOn_MeanXTrial{curr_session}(curr_roi,27:29));
    end
    df_f_aligned_CaOn_MeanXTrialXTime(:,curr_session) = nanmean(df_f_aligned_CaOn_MeanXTrial{curr_session}(:,29:end),2);
    df_f_aligned_CaOn_MeanXTrialXTime_sub(:,curr_session) = nanmean(df_f_aligned_CaOn_MeanXTrial_sub{curr_session}(:,29:end),2);
end

color_map = cbrewer('div','RdBu',64);
color_map = flipud(color_map);
figure; hold on; set(gcf,'color','w','position',[200 200 200*length(df_f_aligned_CaOn) 200])
for curr_session = 1:length(df_f_aligned_CaOn)
    subplot(1,length(df_f_aligned_CaOn),curr_session);
    temp_matrix = df_f_aligned_CaOn_MeanXTrial_sub{curr_session};
    imagesc(temp_matrix,[-1.5 1.5]); colormap(color_map)
    line([29 29],ylim,'color','k','linestyle',':','linewidth',1);
    xticks([1:28:28*5]);
    xticklabels([-1:1:3]);
    xlabel('Time (s)'); ylabel('ROI (axonal segments)');
end
saveas(gcf,'cluster_sub_Heatmap_ActivityAlignedToCaOnset.fig');
saveas(gcf,'cluster_sub_Heatmap_ActivityAlignedToCaOnset.png');

figure; hold on; set(gcf,'color','w','position',[200 200 200 200])
y_limit = [-0.5 1];
for curr_session = 1:length(df_f_aligned_CaOn)
    temp = df_f_aligned_CaOn_MeanXTrial_sub{curr_session};
    temp_11 = nanmean(temp,1);
    temp_12 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1)),1));
    h = area([(temp_11-temp_12)',(2*temp_12)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',repmat([-0.2+curr_session*0.2],1,3),'FaceAlpha',0.3);
    plot(temp_11,'color',repmat([-0.2+curr_session*0.2],1,3),'linewidth',1);
    xlim([1 length(temp_11)]);
    xticks([1:28:28*5]);
    xticklabels([-1:1:3]);
    xlabel('Time (s)'); ylabel('mean df/f');
    ylim(y_limit);
    line([29 29],ylim,'color','k','linestyle',':');
end
axis square;
saveas(gcf,'cluster_sub_Trace_ActivityAlignedToCaOnset.fig');
saveas(gcf,'cluster_sub_Trace_ActivityAlignedToCaOnset.png');

figure; hold on; set(gcf,'color','w','position',[200 200 300 200])
temp = df_f_aligned_CaOn_MeanXTrialXTime_sub;
temp_11 = nanmean(temp,1);
temp_12 = nanstd(temp,[],1)./sqrt(sum(~isnan(temp(:,1)),1));
yyaxis left;
plot(temp_11,'linewidth',1);
for ii = 1:length(df_f_aligned_CaOn)
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
for ii = 1:length(df_f_aligned_CaOn)
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'linewidth',1);
end
ylabel('Norm. mean df/f'); ylim([0 1])
axis square;
saveas(gcf,'cluster_sub_Mean_ActivityAlignedToCaOnset.fig');
saveas(gcf,'cluster_sub_Mean_ActivityAlignedToCaOnset.png');

aaa = df_f_aligned_CaOn_MeanXTrialXTime_sub;
pvalue(1) = signrank(aaa(:,1),aaa(:,2));
pvalue(2) = signrank(aaa(:,1),aaa(:,3));
pvalue(3) = signrank(aaa(:,2),aaa(:,3));
FDR = mafdr(pvalue,'BHFDR', true);

save(['CR_4333447-R_' Date '_jointROI_Analyses.mat'], '-v7.3');
