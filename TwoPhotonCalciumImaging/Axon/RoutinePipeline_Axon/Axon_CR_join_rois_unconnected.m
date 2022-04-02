%%
clear
close all

Animal = 'CR_3575265-LR';
Date = '190308';

General_path = fullfile('F:\Data\MotionCorrection\',[Animal '_warp'],Date,'Rec_1\Z1\motioncorrected_tiff');
load([General_path filesep 'suite2p\plane0\Fall.mat']);
is_cell_index = logical(iscell(:,1));
stat = stat(is_cell_index);
image_template = ones(ops.Ly,ops.Lx);
% get rid of edge artifact
image_template([1:10,ops.Ly-9:ops.Ly],:) = 0;
image_template(:,[1:10,ops.Lx-9:ops.Lx]) = 0;

for ii_cell = 1:length(stat)
    temp_roi = zeros(ops.Ly*ops.Lx,1);
    stat{ii_cell}.ipix = stat{ii_cell}.ipix(stat{ii_cell}.ipix>0);
    temp_roi(stat{ii_cell}.ipix) = 1;
    temp_roi = reshape(temp_roi,ops.Lx,ops.Ly);
    temp_roi = temp_roi';
    temp_roi = temp_roi.*image_template;
    polygon.ROI_mask{ii_cell} = logical(temp_roi);
end

% Get rid of noisy part
c_lim = [0.015 0.1];
for ii_cell = 1:length(stat)
    temp_roi = polygon.ROI_mask{ii_cell};
    figure; clear gca;
    ax1 = axes;
    imagesc(ops.max_proj/max(ops.max_proj(:)),c_lim); axis equal;
    set(ax1,'xlim',[1 ops.Lx],'ylim',[1 ops.Ly]);
    ax2 = axes;
    h = imagesc(temp_roi*c_lim(2),c_lim); axis equal;
    alpha_data = 0.2*double(temp_roi);
    set(ax2,'xlim',[1 ops.Lx],'ylim',[1 ops.Ly]);set(h,'AlphaData',alpha_data);
    title([num2str(ii_cell) '/' num2str(length(stat))]);
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'winter')
    title([num2str(ii_cell) '/' num2str(length(stat))]);
    BW = roipoly;
    temp_roi = temp_roi.*BW;
    imagesc(temp_roi); axis equal;
    title([num2str(ii_cell) '/' num2str(length(stat))]);
    pause;
    polygon.ROI_mask{ii_cell} = logical(temp_roi);
    close gcf;
end

image_template = zeros(ops.Ly,ops.Lx)-10;
for ii_cell = 1:length(polygon.ROI_mask)
    image_template(polygon.ROI_mask{ii_cell}) = ii_cell;
end
figure; hold on; colormap jet;
subplot(1,2,1);
imagesc(image_template); axis equal;
xlim([1 ops.Lx]); ylim([1 ops.Ly]);
subplot(1,2,2);
imagesc(ops.max_proj); axis equal;
xlim([1 ops.Lx]); ylim([1 ops.Ly]);

polygon.clusters = {};
polygon.clusters{1} = [6,21];
polygon.clusters{2} = [26,15];
polygon.clusters{3} = [13,4,3,17,18,1,16];
polygon.clusters{4} = [19,28];
polygon.clusters{5} = [22,7];

polygon.combined_ROI = {};
for ii_cluster = 1:length(polygon.clusters)
    temp_roi = zeros(ops.Ly,ops.Lx);
    for ii_roi = 1:length(polygon.clusters{ii_cluster})
        temp_roi = temp_roi+polygon.ROI_mask{polygon.clusters{ii_cluster}(ii_roi)};
    end
    polygon.combined_ROI{ii_cluster} = logical(temp_roi);
end

cluster_roi = cell2mat(polygon.clusters);
if length(cluster_roi) ~= length(unique(cluster_roi))
    warning('Double Check!!!');
    return
end
cluster_roi = ismember([1:length(polygon.ROI_mask)],cluster_roi);

polygon.Cluster_mask = polygon.ROI_mask(~cluster_roi);
polygon.Cluster_mask = [polygon.Cluster_mask,polygon.combined_ROI];

clear zero_index
for ii_cell = 1:length(polygon.Cluster_mask)
    zero_index(ii_cell) = sum(sum(polygon.Cluster_mask{ii_cell}))==0;
end
polygon.Cluster_mask = polygon.Cluster_mask(~zero_index);

image_template = zeros(ops.Ly,ops.Lx)-10;
for ii_cell = 1:length(polygon.Cluster_mask)
    image_template(polygon.Cluster_mask{ii_cell}) = ii_cell;
end
figure; hold on; colormap jet;
subplot(1,2,1);
imagesc(image_template); axis equal;
xlim([1 ops.Lx]); ylim([1 ops.Ly]);
subplot(1,2,2);
imagesc(ops.max_proj); axis equal;
xlim([1 ops.Lx]); ylim([1 ops.Ly]);

figure; hold on;
for ii = 1:length(polygon.Cluster_mask)
    subplot(5,7,ii)
    imagesc(polygon.Cluster_mask{ii});
end

save([General_path filesep 'summed' filesep Animal(1:3) Date Animal(3:end) '_roi.mat'],'polygon','-v7.3')
cd([General_path filesep 'summed']);




