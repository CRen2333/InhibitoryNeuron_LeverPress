%% Select file
clear all
close all
clc

Filter = true;
FigureExe = false;
EyeROI = true;
MovExe = true;
HoleFill = true;

IN = 'ChAT';
Initial = 'CR';
 
Animal = '3633170-L';
Date = '190411';
Session = 'Rec_1';
FileFormat = 'mat';

TargetPath = ['Z:\People\Chi\TwoP_IN\PupilFitting' filesep Initial '_' Animal filesep Date filesep Session];

% Load baseline information
if exist([TargetPath filesep Initial '_' Date '_' Animal '_SelectFramePer100.mat'],'file')
    load([TargetPath filesep Initial '_' Date '_' Animal '_SelectFramePer100.mat'],'-mat');
else
    disp('No SelectFramePer100.mat file detected!');
    return
end

if exist([TargetPath filesep Initial '_' Date '_' Animal '_' Session '_Baseline_PupilDia_CR.mat'],'file')
    load([TargetPath filesep Initial '_' Date '_' Animal '_' Session '_Baseline_PupilDia_CR.mat'],'Eye_ROI_Mask','Baseline_Pupil_Dia','-mat');
else
    disp('No baseline file detected!');
    return
end

% Estimate Eye ROI and boundry for kmeans
if exist('Mouse_mov_filtered','var')
    [row,col] = find(Eye_ROI_Mask);
    Eye_ROI_filtered = Mouse_mov_filtered(min(row):max(row),min(col):max(col),:);
    clear Mouse_mov_filtered
elseif exist('SelectFramePer100','var')
    [row,col] = find(Eye_ROI_Mask);
    for ii = 1:size(SelectFramePer100,3)
        temp_image = medfilt2(SelectFramePer100(:,:,ii),[3 3],'symmetric');%2-D median filtering
        Eye_ROI_filtered(:,:,ii) = temp_image(min(row):max(row),min(col):max(col));
        clear temp_image
    end
    clear SelectFramePer100
else
    disp('Error');
end

% Select seed area
temp_image = nanmean(Eye_ROI_filtered,3);
[Seed_ROI,Seed_ROI_Mask]=CM_CROP_MOVIE(temp_image);

ClusterNum = 3;
disp('K-means clustering');

for ii = 1:size(Eye_ROI_filtered,3)
    clear temp_image temp_area level
    temp_frame = Eye_ROI_filtered(:,:,ii);

    % multithreshold = k-means clusitering, faster
    level = multithresh(temp_frame,ClusterNum);% Multilevel image thresholds using Otsu's method
    temp_image = imquantize(temp_frame,level);
    % Area control
    for kk = ClusterNum+1:-1:1 % start from brightest
        if nanmean(temp_frame(temp_image(:) == kk)) < level(round(length(level)/2)-1) % the selected area must be bright enough
            continue
        end
        if sum((temp_image(:) == kk).*Seed_ROI_Mask(:)) < 500 % the selected area must colocalize with seed area
            continue
        end
        temp_area = sum(temp_image(:) == kk);
        if temp_area > 2000 % the area must be big enough to be used
            Eye_ROI_kmeans(:,:,ii) = temp_image == kk; % create the mask of interested region and store the ROI in Eye_ROI_kmeans
            clear temp_area
            break
        end
    end
end

% plot the summary of selected area to check whether the ROI is proper
figure; hold on;
for ii = 0:10:10*(floor(size(Eye_ROI_filtered,3)/10)-1)
    subplot(ceil((floor(size(Eye_ROI_filtered,3)/10)-1)/5)+1,5,ii/10+1);
    imagesc(Eye_ROI_kmeans(:,:,ii+1),[0 1]);
    colormap(gray); axis equal; axis off;
end

%%
% Set col and row range for all frames
row_num = size(Eye_ROI_kmeans,1);
col_num = size(Eye_ROI_kmeans,2);
frame_num = size(Eye_ROI_kmeans,3);

% preset ROI
% Horizontal
horiz_cross=squeeze(sum(double(Eye_ROI_kmeans),1));
figure; imagesc(horiz_cross); colormap jet;
[~,col_range] = ginput(2);
col_range_max = round(max(col_range));
col_range_min = round(min(col_range));

% Vertical
vert_cross=squeeze(sum(double(Eye_ROI_kmeans),2));
figure; imagesc(vert_cross); colormap pink;
[~,row_range] = ginput(2);
row_range_max = round(max(row_range));
row_range_min = round(min(row_range));

close all

Eye_ROI_kmeans([1:row_range_min,row_range_max:end],:,:) = false;
Eye_ROI_kmeans(:,[1:col_range_min,col_range_max:end],:) = false;

% Refine k-means
for ii = 1:size(Eye_ROI_filtered,3)
    clear row col Area
    if sum(Eye_ROI_kmeans(:,:,ii)) == 0
        continue
    end
    [row,col] = find(Eye_ROI_kmeans(:,:,ii));
    temp_center = [nanmean(row),nanmean(col)];
    for kk = 1:length(row)
        temp_dist = sqrt((row(kk)-temp_center(1))^2+(col(kk)-temp_center(2))^2);
%         if temp_dist > Baseline_Pupil_Dia
        if temp_dist > (col_range_max-col_range_min)/1.8
            Eye_ROI_kmeans(row(kk),col(kk),ii) = false;
        end
    end
    % Fill the hole
    Eye_ROI_kmeans(:,:,ii) = imfill(Eye_ROI_kmeans(:,:,ii),'holes');

    % Conserve, get rid of eyelid shadow
    [L,n] =  bwlabel(Eye_ROI_kmeans(:,:,ii),4);
    for kk = 1:n
        Area(1,kk) = sum(L(:) == kk);
    end
    [~,index] = max(Area);
    Eye_ROI_kmeans(:,:,ii) = L == index;
end

% Dilate to make the edge smoother
se = strel('disk',2);
for ii = 1:size(Eye_ROI_filtered,3)
    Eye_ROI_kmeans(:,:,ii) = imdilate(Eye_ROI_kmeans(:,:,ii),se);
end

break_length = 0;
disp('Edge detecting');
tic
for ii = 1:size(Eye_ROI_kmeans,3)
    if sum(sum(Eye_ROI_kmeans(:,:,ii))) == 0
        Eye_ROI_Edges(:,:,ii) = false(size(Eye_ROI_kmeans(:,:,1)));
        continue
    end
    [row,col] = find(Eye_ROI_kmeans(:,:,ii));
    temp_center = round([nanmean(row),nanmean(col)]);
    Eye_ROI_Edges(:,:,ii) = edge(Eye_ROI_kmeans(:,:,ii),'canny');
    [row,col] = find(Eye_ROI_Edges(:,:,ii));
    D = pdist([row,col],'euclidean');
    minD = round(prctile(D,15));
    Eye_ROI_Edges(temp_center(1)-minD:temp_center(1)+minD,temp_center(2)-minD:temp_center(2)+minD,ii) = false;
%     Eye_ROI_Edges(min(row_RFL):max(row_RFL),min(col_RFL):max(col_RFL),ii) = false;
    
    % Get rid of relection
    [Gmag,Gdir] = imgradient(Eye_ROI_kmeans(:,:,ii),'prewitt');
    x = (1:size(Eye_ROI_kmeans(:,:,ii),2)) - temp_center(2);
    y = (1:size(Eye_ROI_kmeans(:,:,ii),1)) - temp_center(1);
    [X,Y] = meshgrid(x,y);
    angleMatrix = atan2d(-Y, X);
    angleDifference = mod(Gdir-angleMatrix + 180, 360) - 180;
    angleDifference = abs(angleDifference);
    Eye_ROI_Edges(:,:,ii) = Eye_ROI_Edges(:,:,ii).*(angleDifference>140);
    Eye_ROI_Edges(:,temp_center(2)-break_length:temp_center(2)+break_length,ii) = false;
end
toc

% Get rid of unwanted part
sum_edge_ave = logical(sum(double(Eye_ROI_Edges),3));
figure
h1 = imagesc(Eye_ROI_filtered(:,:,end)/255,[0,1]); hold on; colormap gray
h2 = imagesc(~sum_edge_ave);
set(h2,'AlphaData',sum_edge_ave);
BW_up = roipoly; % draw a polygon ROI to exclude the upper eyelid
BW_low = roipoly; % draw a polygon ROI to exclude the lower eyelid
BW_middle = roipoly;
Unwanted_BW = BW_up + BW_low + BW_middle;
Unwanted_BW = logical(Unwanted_BW);

for ii = 1:size(Eye_ROI_Edges,3)
    temp = Eye_ROI_Edges(:,:,ii);
    temp(Unwanted_BW) = false;
	Eye_ROI_Edges(:,:,ii) = reshape(temp,size(Unwanted_BW));
end   

% Hough transform to fit
absMin = 70; % defult = 70, You may need to change this pretty often;
absMax = 400;
params.rotation = 0;
params.rotationSpan = 90;
params.minAspectRatio = 0.85;
params.randomize = 2;

disp('Fitting...');
tic
for ii = 1:size(Eye_ROI_Edges,3)
    [row,col] = find(Eye_ROI_Edges(:,:,ii));
    D = pdist([row,col],'euclidean');
    params.minMajorAxis = max(absMin,round(prctile(D,85)));
    params.maxMajorAxis = min(absMax,round((col_range_max-col_range_min)/0.9));
    bestFits{1,ii} = ellipseDetection(Eye_ROI_Edges(:,:,ii), params);
end
toc

figure;
for ii = 1:size(Eye_ROI_filtered,3)
    clf;
    imagesc(Eye_ROI_filtered(:,:,ii)/255,[0 1]); colormap gray; axis equal; axis off;
    hold on;
    h = imagesc(Eye_ROI_Edges(:,:,ii)); colormap gray; axis equal; axis off;
    set(h,'AlphaData',Eye_ROI_Edges(:,:,ii));
    ellipse(bestFits{1,ii}(1,3),bestFits{1,ii}(1,4),bestFits{1,ii}(1,5)*pi/180,bestFits{1,ii}(1,1),bestFits{1,ii}(1,2),'r');
    ellipse(bestFits{1,ii}(2,3),bestFits{1,ii}(2,4),bestFits{1,ii}(2,5)*pi/180,bestFits{1,ii}(2,1),bestFits{1,ii}(2,2),'g');
    ellipse(bestFits{1,ii}(3,3),bestFits{1,ii}(3,4),bestFits{1,ii}(3,5)*pi/180,bestFits{1,ii}(3,1),bestFits{1,ii}(3,2),'b');
    title([num2str(ii) '/' num2str(size(Eye_ROI_filtered,3))]);
    pause;
end

close all;
cd(TargetPath);
save([Initial '_' Date '_' Animal '_' Session '_Fitting_PreSet'],'Animal','IN','Initial','Date','TargetPath','Unwanted_BW','Seed_ROI_Mask',...
        'Eye_ROI_filtered','bestFits','ClusterNum','break_length','params','absMin','absMax','col_range_max','col_range_min','row_range_max','row_range_min','-v7.3');

% check blink
aaa = nanmedian(Eye_ROI_filtered,3);
for ii = 1:size(Eye_ROI_filtered,3)
    bbb(ii) = corr2(Eye_ROI_filtered(:,:,ii),aaa);
end
figure;plot(bbb)

prompt = {'Enter coorelation threshold (default: 0.7):'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0.7'};
corrThreshold = inputdlg('Enter coorelation threshold (default: 0.7):','Input',1,{'0.7'});
corrThreshold = str2num(corrThreshold{:});

save([Initial '_' Date '_' Animal '_' Session '_Fitting_PreSet'],'corrThreshold','-append');

close all
clear Eye_ROI_filtered Eye_ROI_kmeans Eye_ROI_Edges bestFits bestFits_refine PupilDia_refine


