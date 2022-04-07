%% Plot Figure S2A, masks of common modules from Thy1-GC6s
load('FigureS2A.mat');
figure; set(gcf,'color','w','position',[200 200 200 200]); hold on;
for ii = 1:length(AverageMode)
    imcontour(AverageMode_threshold{ii},1);
end
xlim([0.5,128.5]);ylim([0.5,128.5]);axis off; axis square;
set(gca,'YDir','reverse')

%% Plot Figure S2B, VIP GC6f example images aligned to movement onset
% Example: CR_3183970-L, session 1
load('FigureS2B.mat');
figure; hold on; set(gcf,'pos',[1400,100,800,600]);
selected_frames = [1:3:76];
for ii = 1:26
    temp_frame = ImageAve(:,selected_frames(ii))-nanmean(ImageAve(:,[5:9]),2);
    temp_frame(1:260) = 0;        
    temp_frame = reshape(temp_frame,[128 128]);
    temp_frame = circshift(temp_frame,[0,2]);
    trans_mask = temp_frame==0;
    trans_mask = reshape(trans_mask,[128 128]);
    subplot(5,6,ii);
    h = imagesc(temp_frame,[-0.005 0.035]); colormap jet;
    set(h,'AlphaData',~trans_mask);
    axis square;
    xlim([0.5 128.5]);ylim([0.5 128.5]);
    xticks([]); yticks([]);
end

%% Plot Figure S2C, SOM GC6f example images aligned to movement onset
% Example: CR_3438521-LR, session 4
load('FigureS2C.mat');
figure; hold on; set(gcf,'pos',[1400,100,800,600]);
selected_frames = [1:3:76];
for ii = 1:26
    temp_frame = ImageAve(:,selected_frames(ii))-nanmean(ImageAve(:,[5:9]),2);
    temp_frame(1:260) = 0;        
    temp_frame = reshape(temp_frame,[128 128]);
    temp_frame = circshift(temp_frame,[0,-1]);
    trans_mask = temp_frame==0;
    trans_mask = reshape(trans_mask,[128 128]);
    subplot(5,6,ii);
    h = imagesc(temp_frame,[-0.01 0.04]); colormap jet;
    set(h,'AlphaData',~trans_mask);
    axis square;
    xlim([0.5 128.5]);ylim([0.5 128.5]);
    xticks([]); yticks([]);
end

%% Plot Figure S2D, PV GC6f example images aligned to movement onset
% Example: CR_3233232-O, session 1
load('FigureS2D.mat');
figure; hold on; set(gcf,'pos',[1400,100,800,600]);
selected_frames = [1:3:76];
for ii = 1:26
    temp_frame = ImageAve(:,selected_frames(ii))-nanmean(ImageAve(:,[5:9]),2);
    temp_frame(1:260) = 0;        
    temp_frame = reshape(temp_frame,[128 128]);
    temp_frame = circshift(temp_frame,[0,-1]);
    trans_mask = temp_frame==0;
    trans_mask = reshape(trans_mask,[128 128]);
    subplot(5,6,ii);
    h = imagesc(temp_frame,[-0.005 0.03]); colormap jet;
    set(h,'AlphaData',~trans_mask);
    axis square;
    xlim([0.5 128.5]);ylim([0.5 128.5]);
    xticks([]); yticks([]);
end

%% Plot Figure S2E, example ICA modes
% VIP_GC6f: 3438544-R; VIP_GFP: CR_4383182-L
load('FigureS2E.mat')
figure; set(gcf,'color','w');
for ii = 1:5
    subaxis(1,5,ii,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = GC6f.module(:,ii);
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end
figure; set(gcf,'color','w');
for ii = 1:5
    subaxis(1,5,ii,'Spacing',0.02,'Padding',0,'Margin',0.01);
    clims = [-3 10];
    image = GFP.module(:,ii);
    imagesc(reshape(image, [128 128]), clims);
    colormap jet;
    axis square;
    axis off;
end

%% Plot Figure S2F, GFP example images aligned to movement onset
% Example: CR_4383182-L, session 2
load('FigureS2F.mat');
figure; hold on; set(gcf,'pos',[1400,100,800,600]);
selected_frames = [1:3:76];
for ii = 1:26
    temp_frame = ImageAve(:,selected_frames(ii))-nanmean(ImageAve(:,[5:9]),2);
    temp_frame(1:260) = 0;        
    temp_frame = reshape(temp_frame,[128 128]);
    temp_frame = circshift(temp_frame,[-15,0]);
    trans_mask = temp_frame==0;
    trans_mask = reshape(trans_mask,[128 128]);
    subplot(5,6,ii);
    h = imagesc(temp_frame,[-0.005 0.035]); colormap jet;
    set(h,'AlphaData',~trans_mask);
    axis square;
    xlim([0.5 128.5]);ylim([0.5 128.5]);
    xticks([]); yticks([]);
end

%% Plot Figure S2G
load('FigureS2GS2H.mat')
ROI = VIP.ROI;
Ordered_ROI = VIP.Ordered_ROI;
reordered_module = VIP.reordered_module;
figure;
set(gcf,'position',[50,50,600,600]);
for curr_ROI = 1:length(VIP.ROI)
    subplot(4,4,curr_ROI); hold on;
    for ii = 1
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        temp_2 = nanstd(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})./sqrt(sum(~isnan(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii})));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor','b','FaceAlpha',0.3);
    end
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    temp_2 = nanstd(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})./sqrt(sum(~isnan(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)})));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.3);
        
    for ii = 1
        temp_1 = nanmean(VIP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)}{ii});
        plot([1:76],temp_1,'color','b','LineWidth',1);
    end    
    temp_1 = nanmean(VIP_GFP.Cued_subtract_ROI_df_f_averageTrace{reordered_module(curr_ROI)});
    plot([1:76],temp_1,'color',[0.7 0.7 0.7],'LineWidth',1);
    
    title(VIP.Ordered_ROI{curr_ROI},'fontsize',8);
    xlim([1 76]); ylabel('Mean df/f'); xlabel('Time (sec)');
    set(gca,'XTick',[16,46,76],'XTickLabel',{'0','1','2'});
    ylim([-0.015, 0.022]); 
    set(gca,'YTick',[-0.02:0.01:0.02],'YTickLabel',[-0.02:0.01:0.02]);
    line([16 16],ylim,'color','k','LineWidth',1,'LineStyle',':');
    axis square;
end

%% Plot Figure S2H
load('FigureS2GS2H.mat')
ROI = VIP.ROI;
Ordered_ROI = VIP.Ordered_ROI;
reordered_module = VIP.reordered_module;
% norm, average across animals
VIP.PostActivity_corticalMean = nanmean(VIP.PostActivity{1},2);
VIP_GFP.PostActivity_corticalMean = nanmean(VIP_GFP.PostActivity,2);
naive_activity = nanmean(VIP.PostActivity_corticalMean);
figure('position',[100,100,200,200],'Color','w')
hold on;
bar([1],nanmean(VIP.PostActivity_corticalMean/naive_activity),0.7,'FaceColor','b','LineStyle','none');
bar([2],nanmean(VIP_GFP.PostActivity_corticalMean/naive_activity),0.7,'FaceColor',[0.7 0.7 0.7],'LineStyle','none');
temp_1 = nanmean(VIP.PostActivity_corticalMean/naive_activity);
temp_2 = nanstd(VIP.PostActivity_corticalMean/naive_activity)/sum(~isnan(VIP.PostActivity_corticalMean))^0.5;
line([1,1],[temp_1-temp_2, temp_1+temp_2],'Color','b','LineWidth',1);
temp_1 = nanmean(VIP_GFP.PostActivity_corticalMean/naive_activity);
temp_2 = nanstd(VIP_GFP.PostActivity_corticalMean/naive_activity)/(sum(~isnan(VIP_GFP.PostActivity_corticalMean))^0.5);
line([2,2],[temp_1-temp_2, temp_1+temp_2],'Color',[0.7 0.7 0.7],'LineWidth',1);
xlim([0.3,2.7]); ylim([-0 1.2])
set(gca,'FontSize',8,'XTick',[1,2],'XTickLabel',{'GCaMP6f','GFP'});
ylabel('Norm. df/f')
