%% Track pupil
%% Select file
clear all
close all
clc

% parameters
Filter = true;
FigureExe = false;
EyeROI = true;
MovExe = true;
HoleFill = true;
ClusterNum = 7; % Defult

IN = 'SOM';
Initial = 'CR';
Animal = '3526549-L';
Date = '180910';

if ispc
    TargetPath = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Pupil' filesep Date];
elseif isunix
    addpath(genpath('/usr/local/lab/People/Chi/WFLP_IN/WFIN_IN/Code/CR_Pupil_Imaging'));
    TargetPath = ['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'Pupil' filesep Date];
end
cd(TargetPath);

% Load baseline information
if exist([Initial '_' Date '_' Animal '_SelectFramePer100.mat'],'file')
    load([Initial '_' Date '_' Animal '_SelectFramePer100.mat'],'-mat');
else
    disp('No SelectFramePer100.mat file detected!');
    return
end

if exist([Initial '_' Date '_' Animal '_Baseline_PupilDia_CR.mat'],'file')
    load([Initial '_' Date '_' Animal '_Baseline_PupilDia_CR.mat'],'Eye_ROI_Mask','Baseline_Pupil_Dia','-mat');
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
        temp_image = medfilt2(SelectFramePer100(:,:,ii),[3 3],'symmetric');
        Eye_ROI_filtered(:,:,ii) = temp_image(min(row):max(row),min(col):max(col));
        clear temp_image
    end
    clear SelectFramePer100
else
    disp('Error');
end

disp('K-means clustering');
ClusterNum = 10;
for ii = 1:size(Eye_ROI_filtered,3)
    clear temp_image temp_area level
    temp_image = Eye_ROI_filtered(:,:,ii);

    % multithreshold = k-means clusitering, faster
    level = multithresh(temp_image,ClusterNum);
    temp_image = imquantize(temp_image,level);
    % Arear control
    for kk = 1:ClusterNum
        temp_area = sum(temp_image(:) == kk);
        if temp_area > 5000
            Eye_ROI_kmeans(:,:,ii) = temp_image == kk;
            clear temp_area
            break
        end
    end
end

figure; hold on;
for ii = 0:30:30*(floor(size(Eye_ROI_filtered,3)/30)-1)
    subplot(7,5,ii/30+1);
    imagesc(Eye_ROI_kmeans(:,:,ii+1),[0 1]);
    colormap(gray); axis equal; axis off;
end

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

% points_num = inputdlg('Enter points number for horiz (2-40):','Input',[1 40]);
% points_num = str2num(points_num{:});
% 
% upside = [frame_num/points_num:frame_num/points_num:frame_num;ones(1,points_num)]';
% downside = [frame_num:-frame_num/points_num:frame_num/points_num;col_num*ones(1,points_num)]';
% tt = [1 1; upside; frame_num col_num/2; downside; 1 col_num; 1 col_num/2];
% hpoly_hor = impoly(gca,tt);
% wait(hpoly_hor);
% BW_dia_hor = createMask(hpoly_hor);
% BW_dia_hor = double(BW_dia_hor);
% BW_dia_hor(logical(BW_dia_hor)~=1)=nan;
% image_to_diam_hor=horiz_cross.*BW_dia_hor;
% figure; imagesc(image_to_diam_hor);

% Vertical
vert_cross=squeeze(sum(double(Eye_ROI_kmeans),2));
figure; imagesc(vert_cross); colormap pink;
[~,row_range] = ginput(2);
row_range_max = round(max(row_range));
row_range_min = round(min(row_range));

% points_num = inputdlg('Enter points number for vert (2-40):','Input',[1 40]);
% points_num = str2num(points_num{:});
% 
% upside = [frame_num/points_num:frame_num/points_num:frame_num;ones(1,points_num)]';
% downside = [frame_num:-frame_num/points_num:frame_num/points_num;row_num*ones(1,points_num)]';
% tt_vert = [1 1; upside; frame_num row_num/2; downside; 1 row_num; 1 row_num/2];
% hpoly_ver = impoly(gca,tt_vert);
% wait(hpoly_ver);
% BW_dia_ver = createMask(hpoly_ver);
% BW_dia_ver=double(BW_dia_ver);
% BW_dia_ver(logical(BW_dia_ver)~=1)=nan;
% image_to_diam_ver=vert_cross.*BW_dia_ver;
% figure; imagesc(image_to_diam_ver);

close all

% for ii = 1:frame_num
%     col_range(ii,1) = find(BW_dia_hor(:,ii)==1,1,'first');
%     col_range(ii,2) = find(BW_dia_hor(:,ii)==1,1,'last');
%     row_range(ii,1) = find(BW_dia_ver(:,ii)==1,1,'first');
%     row_range(ii,2) = find(BW_dia_ver(:,ii)==1,1,'last');
% end

% Refine k-means
% col_range_max = max(col_range(:,2));
% col_range_min = min(col_range(:,1));
% row_range_max = max(row_range(:,2));
% row_range_min = min(row_range(:,1));
Eye_ROI_kmeans([1:row_range_min,row_range_max:end],:,:) = false;
Eye_ROI_kmeans(:,[1:col_range_min,col_range_max:end],:) = false;

% Refine k-means
for ii = 1:size(Eye_ROI_filtered,3)
    clear row col Area
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

% Get rid of reflection
% [~,RFL_ROI_Mask]=CM_CROP_MOVIE(nanmean(double(Eye_ROI_kmeans),3)); % limit to the eye
% [row_RFL,col_RFL] = find(RFL_ROI_Mask);

break_length = 0;
disp('Edge detecting');
tic
for ii = 1:size(Eye_ROI_kmeans,3)
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
    
%     for kk = 1:length(col)
%         if col(kk) < temp_center(2)
%             if Eye_ROI_kmeans(row(kk),col(kk)-1,ii) == 1 &&  Eye_ROI_kmeans(row(kk),col(kk)+1,ii) == 0
%                 Eye_ROI_Edges(row(kk),col(kk),ii) = false;
%             end
%         elseif col(kk) > temp_center(2)
%             if Eye_ROI_kmeans(row(kk),col(kk)+1,ii) == 1 &&  Eye_ROI_kmeans(row(kk),col(kk)-1,ii) == 0
%                 Eye_ROI_Edges(row(kk),col(kk),ii) = false;
%             end
%         end
%         if row(kk) < temp_center(1)
%             if Eye_ROI_kmeans(row(kk)-1,col(kk),ii) == 1 &&  Eye_ROI_kmeans(row(kk)+1,col(kk),ii) == 0
%                 Eye_ROI_Edges(row(kk),col(kk),ii) = false;
%             end
%         elseif row(kk) > temp_center(1)
%             if Eye_ROI_kmeans(row(kk)+1,col(kk),ii) == 1 &&  Eye_ROI_kmeans(row(kk)-1,col(kk),ii) == 0
%                 Eye_ROI_Edges(row(kk),col(kk),ii) = false;
%             end
%         end
%     end
    % ignore bottom and up edge
    Eye_ROI_Edges(:,temp_center(2)-break_length:temp_center(2)+break_length,ii) = false;
end
toc

% Get rid of unwanted part
sum_edge_ave = logical(sum(double(Eye_ROI_Edges),3));
figure
h1 = imagesc(Eye_ROI_filtered(:,:,1)/255,[0,1]); hold on; colormap gray
h2 = imagesc(sum_edge_ave);
set(h2,'AlphaData',sum_edge_ave);
BW_up = roipoly;
BW_mid = roipoly;
BW_low = roipoly;
Unwanted_BW = BW_up + BW_mid + BW_low;
Unwanted_BW = logical(Unwanted_BW);

for ii = 1:size(Eye_ROI_Edges,3)
    temp = Eye_ROI_Edges(:,:,ii);
    temp(Unwanted_BW) = false;
	Eye_ROI_Edges(:,:,ii) = reshape(temp,size(Unwanted_BW));
end   

% Hough transform to fit
absMin = 70; % defult = 70, You may need to change this pretty often;
params.minMajorAxis = 250;
params.maxMajorAxis = 330;
params.rotation = 0;
params.rotationSpan = 90;
params.minAspectRatio = 0.95;
params.randomize = 2;
disp('Fitting...');
tic
for ii = 1:size(Eye_ROI_Edges,3)
    [row,col] = find(Eye_ROI_Edges(:,:,ii));
    D = pdist([row,col],'euclidean');
    params.minMajorAxis = max(absMin,round(prctile(D,85)));
%     params.maxMajorAxis = min(330,round(prctile(D,100)));
    params.maxMajorAxis = min(330,round((col_range_max-col_range_min)/0.9));
    bestFits{1,ii} = ellipseDetection(Eye_ROI_Edges(:,:,ii), params);
end
toc

figure;
for ii = 1:size(Eye_ROI_Edges,3)
    clf;
    imagesc(Eye_ROI_filtered(:,:,ii)/255,[0 1]); colormap gray; axis equal; axis off;
    hold on;
    h = imagesc(Eye_ROI_Edges(:,:,ii)); colormap gray; axis equal; axis off;
    set(h,'AlphaData',Eye_ROI_Edges(:,:,ii));
    ellipse(bestFits{1,ii}(1,3),bestFits{1,ii}(1,4),bestFits{1,ii}(1,5)*pi/180,bestFits{1,ii}(1,1),bestFits{1,ii}(1,2),'r');
    ellipse(bestFits{1,ii}(2,3),bestFits{1,ii}(2,4),bestFits{1,ii}(2,5)*pi/180,bestFits{1,ii}(2,1),bestFits{1,ii}(2,2),'y');
    ellipse(bestFits{1,ii}(3,3),bestFits{1,ii}(3,4),bestFits{1,ii}(3,5)*pi/180,bestFits{1,ii}(3,1),bestFits{1,ii}(3,2),'w');
    title([num2str(ii) '/' num2str(size(Eye_ROI_Edges,3))]);
    pause;
end
close all;

cd(TargetPath);
save([Initial '_' Date '_' Animal '_Fitting_PreSet'],'Animal','IN','Initial','Date','TargetPath','Unwanted_BW',...
        'Eye_ROI_filtered','ClusterNum','break_length','params','col_range_max','col_range_min','row_range_max','row_range_min',...
        'absMin','-v7.3');

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

save([Initial '_' Date '_' Animal '_Fitting_PreSet'],'corrThreshold','-append');

close all
clear Eye_ROI_filtered Eye_ROI_kmeans Eye_ROI_Edges bestFits bestFits_refine PupilDia_refine

