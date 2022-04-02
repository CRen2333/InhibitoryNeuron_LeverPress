%% Have fun with huge craniotomy
% Apply mask of common modules from Thy1GC6 to IN lines, expand from 18x to 25x

%% Check where those modules are
close all
clear all
clc

cd('Z:\People\Chi\WFLP_IN\Coordinates_FromThy1GC6');
OldPrepBrain = double(imread('03_CR_151027_2790556-O_1(2).tif'));
load('CoordinatesForApplyMask.mat', 'AverageMode')
threshold = 3.5;
for ii = 1:length(AverageMode)
    AverageMode_threshold{ii} = AverageMode{ii} > threshold;
end
AverageMode_threshold_all = AverageMode_threshold{1};
for ii = 2:length(AverageMode)
    AverageMode_threshold_all = AverageMode_threshold_all + AverageMode_threshold{ii};
end

% Plot common mask
figure; set(gcf,'color','w','position',[200 200 200 200]); hold on;
for ii = 1:length(AverageMode)
    imcontour(AverageMode_threshold{ii},1);
end
xlim([0.5,128.5]);ylim([0.5,128.5]);axis off; axis square;
set(gca,'YDir','reverse')

% Plot common mask on Thy1-GC6s brain
OldPrepMask = AverageMode_threshold_all;
OldPrepMask(OldPrepMask==0) = nan;
OldPrepBrain = imresize(OldPrepBrain,0.25);
figure
imagesc(OldPrepBrain/6000,[0 1]); colormap(gray); axis square;
hold on;
h = imagesc(OldPrepMask,[0 1]);
set(h,'AlphaData',0.2);
close all

%% Continue
Initial = 'CR';
IN = 'VIP';
Animal = '3438544-R';

switch Animal
    case '3183959-LL'
        Date = '170528';
    case '3183959-LR'
        Date = '170529';
    case '3183958-L'
        Date = '170620';
    case '3183958-R'
        Date = '170627';
    case '3218181-O'
        Date = '170630';
    case '3218181-L'
        Date = '170630';
    case '3218181-R'
        Date = '170627';
    case '3183935-L'
        Date = '170512';
    case '3183970-L'
        Date = '170519';
    case '3183970-O'
        Date = '170519';
    case '3183972-L'
        Date = '170505';
    case '3184011-LL'
        Date = '170522';
    case '3184010-O'
        Date = '170603';
    case '3184012-R'
        Date = '170613';
    case '3218183-O'
        Date = '170616';
    case '3218183-R'
        Date = '170614';
    case '3161016-O'
        Date = '170418';
    case '3161018-R'
        Date = '170411';
    case '3233232-O'
        Date = '170706';
    case '3233232-L'
        Date = '170706';
    case '3233232-R'
        Date = '170707';
    case '3233233-O'
        Date = '170721';
    case '3233233-L'
        Date = '170721';
    case '3233233-R'
        Date = '170721';
    case '3258531-O'
        Date = '171129';
    case '3258531-L'
        Date = '171129';
    case '3358883-O'
        Date = '171214';
    case '3358883-R'
        Date = '171214';
    case '3358884-L'
        Date = '171213';
    case '3358884-R'
        Date = '171213';
    case '3271320-L'
        Date = '171105';
    case '3333408-O'
        Date = '171105';
    case '3333408-L'
        Date = '171105';
    case '3333408-R'
        Date = '171105';
    case '3373693-O'
        Date = '180111';
    case '3373693-L'
        Date = '180111';
    case '3373693-R'
        Date = '180110';
    case '3373693-LR'
        Date = '180110';
    case '3420509-L'
        Date = '180129';
    case '3420509-R'
        Date = '180129';
    case '3373785-O'
        Date = '180327';
    case '3373785-L'
        Date = '180327';
    case '3373785-R'
        Date = '180327';
    case '3438483-L'
        Date = '180501';
    case '3491479-L'
        Date = '180501';
    case '3491479-LR'
        Date = '180503';
    case '3491479-R'
        Date = '180508';
    case '3438500-O'
        Date = '180508';
    case '3438500-L'
        Date = '180508';
    case '3438544-O'
        Date = '180624';
    case '3438544-L'
        Date = '180624';
    case '3438544-R'
        Date = '180624';
    case '3453262-O'
        Date = '180625';
    case '3453262-L'
        Date = '180625';
    case '3557180-R'
        Date = '180804';
    case '3453312-O'
        Date = '180719';
    case '3453312-L'
        Date = '180719';
    case '3547207-LR'
        Date = '181103';
    case '3575195-O'
        Date = '190226';
    case '3547274-O'
        Date = '190308';
    case '3575194-O'
        Date = '190318';
    case '3575194-L'
        Date = '190319';
    case '4383182-O'
        Date = '210904';
    case '4383182-L'
        Date = '210904';
    case '4383183-O'
        Date = '210904';
end

cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
load([Initial '_' Animal '_RecICA_80.mat'],'sortMode_Retained');
threshold = 3.5;
for curr_ROI = 1:length(sortMode_Retained)
    sortMode_mask{curr_ROI} = sortMode_Retained{curr_ROI}>threshold;
end

cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'df_f\' Date]);

NewPrepBrain = double(imread(['CR_' Date '_' Animal '_01(2).tif']));
NewPrepBrain = imresize(NewPrepBrain,0.25);
load(['CR_' Date '_' Animal '_01(2).pixel'],'-mat');
temp_mask = true(128,128);
temp_mask(roiPixelNum) = false;
NewPrepBrain(temp_mask) = 0;
load(['CR_' Date '_' Animal '_01(2).coordinatePixel'],'-mat');
NewPrepBrain2 = circshift(NewPrepBrain,[71-round(coordinate{1,1}(2)/4),64-round(coordinate{1,1}(1)/4)]);
% NewPrepBrain2 = circshift(NewPrepBrain,[44-round(coordinate{1,1}(2)/4),64-round(coordinate{1,1}(1)/4)]);

NewPrepMask = zeros(128,128);
for curr_ROI = 1:length(sortMode_Retained)
    NewPrepMask = NewPrepMask + sortMode_mask{curr_ROI}*rand(1);
end

figure
h1 = imagesc(NewPrepBrain2/5000,[-0.5 1.5]); colormap(gray); axis square; axis off
transMask = ones(128,128);
transMask(NewPrepBrain2==0) = 0;
set(h1,'AlphaData',transMask);
hold on;
h2 = imagesc(NewPrepMask,[-0.5 1.5]);
transMask = ones(128,128);
transMask(NewPrepMask==0) = 0;
transMask(NewPrepMask~=0) = 0.4;
set(h2,'AlphaData',transMask);

close all
% Expand old prep by 1.33x
OldPrepMask(isnan(OldPrepMask)) = 0;
OldPrepMask_25x = imresize(OldPrepMask,1.33);
OldPrepMask_25x(isnan(OldPrepMask_25x)) = nan;
OldPrepBrain_25x = imresize(OldPrepBrain,1.32);

figure(1)
h = imagesc(OldPrepBrain_25x/6000,[0 1]); colormap(gray); axis square;
% hold on;
% h = imagesc(OldPrepMask_25x,[0 1]);
% set(h,'AlphaData',0.2);
hold on;
imcontour(OldPrepMask_25x,1);

% Extend new prep
NewPrepBrain2_Ext = zeros(size(OldPrepMask_25x));
boundry = floor((size(NewPrepBrain2_Ext,1)-size(NewPrepBrain2,1))./2);
NewPrepBrain2_Ext(boundry:boundry+127,boundry:boundry+127) = NewPrepBrain2;
NewPrepMask_Ext(boundry:boundry+127,boundry:boundry+127) = NewPrepMask;

figure(2)
h1 = imagesc(NewPrepBrain2_Ext/6000,[-0.5 1.5]); colormap(gray); axis square; axis off
transMask = ones(size(NewPrepBrain2_Ext));
transMask(NewPrepBrain2_Ext==0) = 0;
set(h1,'AlphaData',transMask);
% hold on;
% h2 = imagesc(NewPrepMask_Ext,[-0.5 1.5]);
% transMask = ones(size(NewPrepMask_Ext));
% transMask(NewPrepMask_Ext==0) = 0;
% transMask(NewPrepMask_Ext~=0) = 0.4;
% set(h2,'AlphaData',transMask);
hold on;
imcontour(NewPrepMask_Ext,10);

% Overlay the new and old map
figure(9); set(gcf,'position',[100,100,800,800]);
h1 = imagesc(OldPrepBrain_25x/6000,[0 1]); colormap(gray); axis square;
hold on
h2 = imagesc(NewPrepBrain2_Ext/6000,[0 1]); colormap(gray); axis square; axis off
transMask = ones(size(NewPrepBrain2_Ext));
transMask(NewPrepBrain2_Ext==0) = 0;
transMask(NewPrepBrain2_Ext~=0) = 0.4;
set(h2,'AlphaData',transMask);
pause; % don't know why if directly continue the h2 will disappear
hold on;
imcontour(OldPrepMask_25x,1);
pause;
hold on;
imcontour(NewPrepMask_Ext,4);

% Shift
switch Animal
    case '3183959-LL'
        Shift_coordinate = [0,8];
    case '3183959-LR'
        Shift_coordinate = [0,8];
    case '3183958-L'
        Shift_coordinate = [0,8];
    case '3183958-R'
        Shift_coordinate = [0,2];
    case '3218181-O'
        Shift_coordinate = [0,4];
    case '3218181-L'
        Shift_coordinate = [0,8];
    case '3218181-R'
        Shift_coordinate = [0,9];
    case '3183935-L'
        Shift_coordinate = [0,2];
    case '3183970-L'
        Shift_coordinate = [0,10];
    case '3183970-O'
        Shift_coordinate = [0,11];
    case '3183972-L'
        Shift_coordinate = [0,10];
    case '3184011-LL'
        Shift_coordinate = [0,7];
    case '3184010-O'
        Shift_coordinate = [-1,2]; % Bregma
    case '3184012-R'
        Shift_coordinate = [0,4];
    case '3218183-O'
        Shift_coordinate = [0,7];
    case '3218183-R'
        Shift_coordinate = [0,8];
    case '3161016-O'
        Shift_coordinate = [0,10];
    case '3161018-R'
        Shift_coordinate = [0,8];
    case '3233232-O'
        Shift_coordinate = [0,4];
    case '3233232-L'
        Shift_coordinate = [0,6];
    case '3233232-R'
        Shift_coordinate = [0,4];
    case '3233233-O'
        Shift_coordinate = [0,4];
    case '3233233-L'
        Shift_coordinate = [0,8];
    case '3233233-R'
        Shift_coordinate = [0,12];
    case '3258531-O'
        Shift_coordinate = [2,6];
    case '3258531-L'
        Shift_coordinate = [1,8];
    case '3358883-O'
        Shift_coordinate = [0,8];
    case '3358883-R'
        Shift_coordinate = [0,7];
    case '3358884-L'
        Shift_coordinate = [0,8];
    case '3358884-R'
        Shift_coordinate = [0,4];
    case '3271320-L'
        Shift_coordinate = [0,6];
    case '3333408-O'
        Shift_coordinate = [0,6];
    case '3333408-L'
        Shift_coordinate = [0,10];
    case '3333408-R'
        Shift_coordinate = [0,7];
    case '3373693-O'
        Shift_coordinate = [0,10];
    case '3373693-L'
        Shift_coordinate = [-2,12];
    case '3373693-R'
        Shift_coordinate = [0,12];
    case '3373693-LR'
        Shift_coordinate = [0,14];
    case '3240509-L'
        Shift_coordinate = [0,10];
    case '3240509-R'
        Shift_coordinate = [0,14];
    case '3373785-O'
        Shift_coordinate = [0,16];
    case '3373785-L'
        Shift_coordinate = [-1,16];
    case '3373785-R'
        Shift_coordinate = [0,8];
    case '3438483-L'
        Shift_coordinate = [-1,8];
    case '3491479-L'
        Shift_coordinate = [0,8];
    case '3491479-LR'
        Shift_coordinate = [0,10];
    case '3491479-R'
        Shift_coordinate = [0,4];
    case '3438500-O'
        Shift_coordinate = [0,8];
    case '3438500-L'
        Shift_coordinate = [0,4];
    case '3438544-O'
        Shift_coordinate = [0,6];
    case '3438544-L'
        Shift_coordinate = [0,10];
    case '3438544-R'
        Shift_coordinate = [0,12];
    case '3453262-O'
        Shift_coordinate = [0,12];
    case '3453262-L'
        Shift_coordinate = [0,12];
    case '3557180-R'
        Shift_coordinate = [0,10];
    case '3453312-O'
        Shift_coordinate = [0,16];
    case '3453312-L'
        Shift_coordinate = [0,10];
    case '3547207-LR'
        Shift_coordinate = [0,10];
    case '3575195-O'
        Shift_coordinate = [0,10];
    case '3575194-O'
        Shift_coordinate = [0,10];
    case '3575194-L'
        Shift_coordinate = [0,10];
end
temp_NewPrepBrain2_Ext = circshift(NewPrepBrain2_Ext,[Shift_coordinate(2),Shift_coordinate(1)]);
temp_NewPrepMask_Ext = circshift(NewPrepMask_Ext,[Shift_coordinate(2),Shift_coordinate(1)]);

figure(4); set(gcf,'position',[100,100,800,800]);
h1 = imagesc(OldPrepBrain_25x/5000,[0 1]); colormap(gray); axis square;
hold on
h2 = imagesc(temp_NewPrepBrain2_Ext/6000,[0 1]); colormap(gray); axis square; axis off
transMask = ones(size(temp_NewPrepBrain2_Ext));
transMask(temp_NewPrepBrain2_Ext==0) = 0;
transMask(temp_NewPrepBrain2_Ext~=0) = 0.5;
set(h2,'AlphaData',transMask);

figure(5)
h1 = imagesc(OldPrepBrain_25x/6000,[0 1]); colormap(gray); axis square;
hold on;
imcontour(OldPrepMask_25x,1);
hold on;
% imcontour(temp_NewPrepMask_Ext,10);
h2 = imagesc(temp_NewPrepMask_Ext,[-0.5 1.5]);
transMask = ones(size(temp_NewPrepMask_Ext));
transMask(temp_NewPrepMask_Ext==0) = 0;
transMask(temp_NewPrepMask_Ext~=0) = 0.4;
set(h2,'AlphaData',transMask);

% Shift old mask
switch Animal
    case '3183959-LL'
        Shift_coordinate2 = [-1,-8];
    case '3183959-LR'
        Shift_coordinate2 = [0,-6];
    case '3183958-L'
        Shift_coordinate2 = [-1,-8];
    case '3183958-R'
        Shift_coordinate2 = [-1,-4];
    case '3218181-O'
        Shift_coordinate2 = [-1,-2];
    case '3218181-L'
        Shift_coordinate2 = [-5,-8];
    case '3218181-R'
        Shift_coordinate2 = [-2,-9];
    case '3183935-L'
        Shift_coordinate2 = [-2,-2];
    case '3183970-L'
        Shift_coordinate2 = [-1,-10];
    case '3183970-O'
        Shift_coordinate2 = [0,-11];
    case '3183972-L'
        Shift_coordinate2 = [0,-10];
    case '3184011-LL'
        Shift_coordinate2 = [-1,-7];
    case '3184010-O'
        Shift_coordinate2 = [0,-2];
    case '3184012-R'
        Shift_coordinate2 = [-3,-4];
    case '3218183-O'
        Shift_coordinate2 = [-2,-7];
    case '3218183-R'
        Shift_coordinate2 = [-5,-8];
    case '3161016-O'
        Shift_coordinate2 = [-3,-10];
    case '3161018-R'
        Shift_coordinate2 = [-1,-8];
    case '3233232-O'
        Shift_coordinate2 = [1,-4];
    case '3233232-L'
        Shift_coordinate2 = [-1,-6];
    case '3233232-R'
        Shift_coordinate2 = [-1,-4];
    case '3233233-O'
        Shift_coordinate2 = [2,-4];
    case '3233233-L'
        Shift_coordinate2 = [2,-8];
    case '3233233-R'
        Shift_coordinate2 = [2,-12];
    case '3258531-O'
        Shift_coordinate2 = [-2,-6];
    case '3258531-L'
        Shift_coordinate2 = [-2,-7];
    case '3358883-O'
        Shift_coordinate2 = [-1,-7];
    case '3358883-R'
        Shift_coordinate2 = [1,-6];
    case '3358884-L'
        Shift_coordinate2 = [-2,-10];
    case '3358884-R'
        Shift_coordinate2 = [-1,-8];
    case '3271320-L'
        Shift_coordinate2 = [-1,-7];
    case '3333408-O'
        Shift_coordinate2 = [-1,-6];
    case '3333408-L'
        Shift_coordinate2 = [-5,-10];
    case '3333408-R'
        Shift_coordinate2 = [-2,-6];
    case '3373693-O'
        Shift_coordinate2 = [-1,-10];
    case '3373693-L'
        Shift_coordinate2 = [-1,-12];
    case '3373693-R'
        Shift_coordinate2 = [-1,-12];
    case '3373693-LR'
        Shift_coordinate2 = [-1,-14];
    case '3240509-L'
        Shift_coordinate2 = [-1,-10];
    case '3240509-R'
        Shift_coordinate2 = [-1,-14];
    case '3373785-O'
        Shift_coordinate2 = [-1,-14];
    case '3373785-L'
        Shift_coordinate2 = [1,-14];
    case '3373785-R'
        Shift_coordinate2 = [-1,-10];
    case '3438483-L'
        Shift_coordinate2 = [1,-6];
    case '3491479-L'
        Shift_coordinate2 = [-2,-10];
    case '3491479-LR'
        Shift_coordinate2 = [1,-10];
    case '3491479-R'
        Shift_coordinate2 = [-3,-6];
    case '3438500-O'
        Shift_coordinate2 = [-1,-8];
    case '3438500-L'
        Shift_coordinate2 = [-2,-4];
    case '3438544-O'
        Shift_coordinate2 = [-2,-10];
    case '3438544-L'
        Shift_coordinate2 = [-1,-10];
    case '3438544-R'
        Shift_coordinate2 = [-2,-14];
    case '3453262-O'
        Shift_coordinate2 = [-2,-12];
    case '3453262-L'
        Shift_coordinate2 = [-1,-10];
    case '3557180-R'
        Shift_coordinate2 = [1,-11];
    case '3453312-O'
        Shift_coordinate2 = [1,-14];
    case '3453312-L'
        Shift_coordinate2 = [-3,-14];
    case '3547207-LR'
        Shift_coordinate2 = [-2,-6];
    case '3575195-O'
        Shift_coordinate2 = [0,-6];
    case '3575194-O'
        Shift_coordinate2 = [-1,-10];
    case '3575194-L'
        Shift_coordinate2 = [-1,-8];
    case '4383182-O'
        Shift_coordinate2 = [-1,3];
    case '4383182-L'
        Shift_coordinate2 = [-2,9];
    case '4383183-O'
        Shift_coordinate2 = [-2,9];
end 
OldPrepMask_25x_shifted = circshift(OldPrepMask_25x,[Shift_coordinate2(2),Shift_coordinate2(1)]);
figure(6)
h1 = imagesc(NewPrepBrain2_Ext/9000,[0 1]); colormap(gray); axis square;
hold on;
imcontour(OldPrepMask_25x_shifted,1);
% h2 = imagesc(OldPrepMask_25x_shifted); colormap(gray); axis square; axis off
% transMask = ones(size(OldPrepMask_25x_shifted));
% transMask(OldPrepMask_25x_shifted==0) = 0;
% set(h2,'AlphaData',transMask);

% Get mask for blood vessel
switch Animal
    case '3183959-LL'
        Threshold = 300;
    case '3183959-LR'
        Threshold = 500;
    case '3183958-L'
        Threshold = 300;
    case '3183958-R'
        Threshold = 300;
    case '3218181-O'
        Threshold = 250;
    case '3218181-L'
        Threshold = 250;
    case '3218181-R'
        Threshold = 250;
    case '3183935-L'
        Threshold = 300;
    case '3183970-L'
        Threshold = 400;
    case '3183970-O'
        Threshold = 250;
    case '3183972-L'
        Threshold = 200;
    case '3184011-LL'
        Threshold = 450;
    case '3184010-O'
        Threshold = 300;
    case '3184012-R'
        Threshold = 300;
    case '3218183-O'
        Threshold = 350;
    case '3218183-R'
        Threshold = 300;
    case '3161016-O'
        Threshold = 200;
    case '3161018-R'
        Threshold = 200;
    case '3233232-O'
        Threshold = 250;
    case '3233232-L'
        Threshold = 300;
    case '3233232-R'
        Threshold = 300;
    case '3233233-O'
        Threshold = 300;
    case '3233233-L'
        Threshold = 250;
    case '3233233-R'
        Threshold = 250;
    case '3258531-O'
        Threshold = 250;
    case '3258531-L'
        Threshold = 250;
    case '3358883-O'
        Threshold = 300;
    case '3358883-R'
        Threshold = 300;
    case '3358884-L'
        Threshold = 300;
    case '3358884-R'
        Threshold = 250;
    case '3271320-L'
        Threshold = 250;
    case '3333408-O'
        Threshold = 250;
    case '3333408-L'
        Threshold = 250;
    case '3333408-R'
        Threshold = 250;
    case '3373693-O'
        Threshold = 250;
    case '3373693-L'
        Threshold = 250;
    case '3373693-R'
        Threshold = 250;
    case '3373693-LR'
        Threshold = 250;
    case '3420509-L'
        Threshold = 500;
    case '3420509-R'
        Threshold = 500;
    case '3373785-O'
        Threshold = 250;
    case '3373785-L'
        Threshold = 250;
    case '3373785-R'
        Threshold = 250;
    case '3438483-L'
        Threshold = 250;
    case '3491479-L'
        Threshold = 250;
    case '3491479-LR'
        Threshold = 250;
    case '3491479-R'
        Threshold = 250;
    case '3438500-O'
        Threshold = 250;
    case '3438500-L'
        Threshold = 250;
    case '3438544-O'
        Threshold = 250;
    case '3438544-L'
        Threshold = 250;
    case '3438544-R'
        Threshold = 250;
    case '3453262-O'
        Threshold = 250;
    case '3453262-L'
        Threshold = 250;
    case '3557180-R'
        Threshold = 300;
    case '3453312-O'
        Threshold = 250;
    case '3453312-L'
        Threshold = 250;
    case '3547207-LR'
        Threshold = 250;
    case '3575195-O'
        Threshold = 250;
    case '3575194-O'
        Threshold = 250;
    case '3575194-L'
        Threshold = 250;
    case '4383182-O'
        Threshold = 250;
    case '4383182-L'
        Threshold = 250;
    case '4383183-O'
        Threshold = 250;
end
bloodVessels = VesselExtract(NewPrepBrain2, Threshold);
figure;
subplot(1,2,1);imagesc(NewPrepBrain2);title('Input Image'); colormap(gray); axis square;
subplot(1,2,2);imagesc(bloodVessels);title('Extracted Blood Vessels'); colormap(gray); axis square;

switch Animal
    case '3183959-LL'
        Threshold2 = 300;
    case '3183959-LR'
        Threshold2 = 500;
    case '3183958-L'
        Threshold2 = 300;
    case '3183958-R'
        Threshold2 = 300;
    case '3218181-O'
        Threshold2 = 250;
    case '3218181-L'
        Threshold2 = 250;
    case '3218181-R'
        Threshold2 = 250;
    case '3183935-L'
        Threshold2 = 200;
    case '3183970-L'
        Threshold2 = 500;
    case '3183970-O'
        Threshold2 = 400;
    case '3183972-L'
        Threshold2 = 250;
    case '3184011-LL'
        Threshold2 = 450;
    case '3184010-O'
        Threshold2 = 300;
    case '3184012-R'
        Threshold2 = 300;
    case '3218183-O'
        Threshold2 = 350;
    case '3218183-R'
        Threshold2 = 300;
    case '3161016-O'
        Threshold2 = 200;
    case '3161018-R'
        Threshold2 = 200;
    case '3233232-O'
        Threshold2 = 250;
    case '3233232-L'
        Threshold2 = 300;
    case '3233232-R'
        Threshold2 = 300;
    case '3233233-O'
        Threshold2 = 300;
    case '3233233-L'
        Threshold2 = 250;
    case '3233233-R'
        Threshold2 = 250;
    case '3258531-O'
        Threshold2 = 250;
    case '3258531-L'
        Threshold2 = 300;
    case '3358883-O'
        Threshold2 = 250;
    case '3358883-R'
        Threshold2 = 250;
    case '3358884-L'
        Threshold2 = 250;
    case '3358884-R'
        Threshold2 = 250;
    case '3271320-L'
        Threshold2 = 250;
    case '3333408-O'
        Threshold2 = 250;
    case '3333408-L'
        Threshold2 = 250;
    case '3333408-R'
        Threshold2 = 250;
    case '3373693-O'
        Threshold2 = 250;
    case '3373693-L'
        Threshold2 = 250;
    case '3373693-R'
        Threshold2 = 250;
    case '3373693-LR'
        Threshold2 = 250;
    case '3420509-L'
        Threshold2 = 500;
    case '3420509-R'
        Threshold2 = 500;
    case '3373785-O'
        Threshold2 = 300;
    case '3373785-L'
        Threshold2 = 250;
    case '3373785-R'
        Threshold2 = 250;
    case '3438483-L'
        Threshold2 = 300;
    case '3491479-L'
        Threshold2 = 250;
    case '3491479-LR'
        Threshold2 = 300;
    case '3491479-R'
        Threshold2 = 300;
    case '3438500-O'
        Threshold2 = 250;
    case '3438500-L'
        Threshold2 = 250;
    case '3438544-O'
        Threshold2 = 250;
    case '3438544-L'
        Threshold2 = 250;
    case '3438544-R'
        Threshold2 = 300;
    case '3453262-O'
        Threshold2 = 300;
    case '3453262-L'
        Threshold2 = 300;
    case '3557180-R'
        Threshold2 = 350;
    case '3453312-O'
        Threshold2 = 300;
    case '3453312-L'
        Threshold2 = 300;
    case '3547207-LR'
        Threshold2 = 250;
    case '3575195-O'
        Threshold2 = 270;
    case '3575194-O'
        Threshold2 = 270;
    case '3575194-L'
        Threshold2 = 300;
    case '4383182-O'
        Threshold2 = 280;
    case '4383182-L'
        Threshold2 = 260;
    case '4383183-O'
        Threshold2 = 260;
end
BloodMask = bloodVessels>Threshold2;
figure
imagesc(BloodMask);

% fill
BloodMask_fill = BloodMask;
for xii = 2:127
    for yii = 2:127
        if BloodMask_fill(xii,yii)
            continue
        end
        temp_square = BloodMask(xii-1:xii+1,yii-1:yii+1);
        temp_sum = sum(temp_square(:))-temp_square(2,2);
        if temp_sum > 4
            BloodMask_fill(xii,yii) = true;
        end
    end
end
figure
imagesc(BloodMask_fill);

% Get rid off dust
dust_value = prctile(NewPrepBrain2(NewPrepBrain2~=0),95);
for xii = 1:128
    for yii = 1:128
        if BloodMask_fill(xii,yii)
            if NewPrepBrain2(xii,yii) > dust_value
                BloodMask_fill(xii,yii) = false;
            end
        end
    end
end
figure
imagesc(BloodMask_fill);

% Shift mask for each ROI
sortMode_Retained{1}(:,1:3) = 0; % get rid of noise
% Boundry_Mask = sortMode_Retained{1} == 0;
Boundry_Mask = NewPrepBrain2 == 0;
for ii = 1:length(AverageMode_threshold)
    temp_mode_ext = imresize(AverageMode_threshold{ii},1.33);
    temp_mode_ext_shift = circshift(temp_mode_ext,[Shift_coordinate2(2),Shift_coordinate2(1)]);
    temp_mode_ext_shift = temp_mode_ext_shift(boundry:boundry+127,boundry:boundry+127);
    temp_modex_ext_shift_BV = temp_mode_ext_shift.*~BloodMask_fill;
    Final_Mask{ii} = logical(temp_modex_ext_shift_BV.*~Boundry_Mask);
    Final_Mask_withBV{ii} = logical(temp_mode_ext_shift.*~Boundry_Mask);
end

Final_Mask_all = zeros(128,128);
for curr_ROI = 1:length(Final_Mask)
    Final_Mask_all = Final_Mask_all + Final_Mask{curr_ROI}*rand(1);
end

Final_Mask_withBV_all = zeros(128,128);
for curr_ROI = 1:length(Final_Mask)
    Final_Mask_withBV_all = Final_Mask_withBV_all + Final_Mask_withBV{curr_ROI}*rand(1);
end

figure
h1 = imagesc(NewPrepBrain2/6000,[-0.5 1.5]); colormap(gray); axis square; axis off
transMask = ones(size(NewPrepBrain2));
transMask(NewPrepBrain2==0) = 0;
set(h1,'AlphaData',transMask);
hold on;
h2 = imagesc(Final_Mask_all,[-0.5 1.5]);
transMask = ones(size(Final_Mask_all));
transMask(Final_Mask_all==0) = 0;
transMask(Final_Mask_all~=0) = 0.8;
set(h2,'AlphaData',transMask);

figure
h1 = imagesc(NewPrepBrain2/6000,[-0.5 1.5]); colormap(gray); axis square; axis off
transMask = ones(size(NewPrepBrain2));
transMask(NewPrepBrain2==0) = 0;
set(h1,'AlphaData',transMask);
hold on;
h2 = imagesc(Final_Mask_withBV_all,[-0.5 1.5]);
transMask = ones(size(Final_Mask_withBV_all));
transMask(Final_Mask_withBV_all==0) = 0;
transMask(Final_Mask_withBV_all~=0) = 0.8;
set(h2,'AlphaData',transMask);

% Exclude overlapping pixels
temp_mask_all = Final_Mask{1};
for roi_ii = 2:length(Final_Mask)
    temp_mask_all = temp_mask_all + Final_Mask{roi_ii};
end
temp_mask_all = temp_mask_all <= 1;
for roi_ii = 1:length(Final_Mask)
    Final_Mask_OP{roi_ii} = logical(Final_Mask{roi_ii}.*temp_mask_all);
end

temp_mask_all = Final_Mask_withBV{1};
for roi_ii = 2:length(Final_Mask)
    temp_mask_all = temp_mask_all + Final_Mask_withBV{roi_ii};
end
temp_mask_all = temp_mask_all <= 1;
for roi_ii = 1:length(Final_Mask)
    Final_Mask_withBV_OP{roi_ii} = logical(Final_Mask_withBV{roi_ii}.*temp_mask_all);
end

temp_mask_all = Final_Mask_OP{1};
for roi_ii = 2:length(Final_Mask_OP)
    temp_mask_all = temp_mask_all + Final_Mask_OP{roi_ii};
end
figure;
imagesc(temp_mask_all); colormap gray; axis square; axis off;
pause
close(gcf);

temp_mask_all = Final_Mask_withBV_OP{1};
for roi_ii = 2:length(Final_Mask_withBV_OP)
    temp_mask_all = temp_mask_all + Final_Mask_withBV_OP{roi_ii};
end
figure;
imagesc(temp_mask_all); colormap gray; axis square; axis off;
pause
close(gcf);

% cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Pharm' filesep 'EventAligned_Gap500']);
cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);

save([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Shift_coordinate2','Final_Mask','Final_Mask_OP','Final_Mask_withBV','Final_Mask_withBV_OP','-v7.3');

%% For animals that already have masks, exclude overelapping pixels if hasn't yet
clear all
close all
clc

Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O'...
     '3271320-L','3333408-L','3333408-O','3333408-R','3373693-O','3373693-L','3373693-R','3373693-LR'};
IN = 'VIP';

% Exclude overlapping pixels
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    load([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Final_Mask','-mat');
    temp_mask_all = Final_Mask{1};
    for roi_ii = 2:length(Final_Mask)
        temp_mask_all = temp_mask_all + Final_Mask{roi_ii};
    end
    temp_mask_all = temp_mask_all <= 1;
    for roi_ii = 1:length(Final_Mask)
        Final_Mask_OP{roi_ii} = logical(Final_Mask{roi_ii}.*temp_mask_all);
    end

    temp_mask_all = Final_Mask_OP{1};
    for roi_ii = 2:length(Final_Mask_OP)
        temp_mask_all = temp_mask_all + Final_Mask_OP{roi_ii};
    end
    figure;
    imagesc(temp_mask_all); colormap gray; axis square; axis off;
    pause
    close(gcf);

    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    save([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Final_Mask_OP','-append');
    clear Final_Mask Final_Mask_OP
end

%% Get df_f for each cortical region
clear all
close all
clc

IN = 'VIP';
Initial = 'CR';
Animals = {'3438544-O','3438544-L','3438544-R'};

OP = true;

for curr_animal = 3:length(Animals)
    clearvars -except IN Initial Animals OP curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500\MovOnsetAligned']);
    All_file = dir(cd);
    All_filename = {All_file.name};
    All_filename = All_filename(3:end)';
    for curr_session = 1:length(All_filename)
        load(All_filename{curr_session},'');
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
    end
    clear AlignedIm_MovOnset
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500']);
    load([Initial '_' Animal '_RecICA_80.mat'],'RecICA_Mov');
    if OP
        load([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Final_Mask_OP','-mat');
        Final_Mask = Final_Mask_OP;
    else
        load([Initial '_' Animal '_Mask_from_Thy1GC6s'],'Final_Mask','-mat');
    end

    % Get df/f for each ROI
    for curr_ROI = 1:length(Final_Mask)
        for curr_session = 1:min(11,length(RecICA_Mov))
            if isempty(RecICA_Mov{curr_session})
                df_f_ROI{curr_ROI}{curr_session} = [];
            else
                df_f_ROI{curr_ROI}{curr_session} = nanmean(RecICA_Mov{curr_session}(Final_Mask{curr_ROI},:),1);
            end
            if ismember(Animal,{'3438500-O','3438544-O','3438544-L'})
                trialnum{curr_session} = size(RecICA_Mov{curr_session},2)/136;
                df_f_ROI_trial_arranged{curr_ROI}{curr_session} = reshape(df_f_ROI{curr_ROI}{curr_session},136,trialnum{curr_session})';
            else
                trialnum{curr_session} = size(RecICA_Mov{curr_session},2)/76;
                df_f_ROI_trial_arranged{curr_ROI}{curr_session} = reshape(df_f_ROI{curr_ROI}{curr_session},76,trialnum{curr_session})';
            end
            index = logical(Index_Cued_trial{curr_session}.*~Index_Catch_trial{curr_session});
            if length(index) ~= trialnum{curr_session}
                disp(['Index:' num2str(length(index)) ', TrailNum: ' num2str(trialnum{curr_session})]);
                index = index(1:trialnum{curr_session});
            end
            df_f_ROI_trial_arranged{curr_ROI}{curr_session} = df_f_ROI_trial_arranged{curr_ROI}{curr_session}(index,1:76);
            df_f_ROI_sub_1st_trial_arranged{curr_ROI}{curr_session} = bsxfun(@minus,df_f_ROI_trial_arranged{curr_ROI}{curr_session},df_f_ROI_trial_arranged{curr_ROI}{curr_session}(:,1));
            df_f_ROI_sub_16th_trial_arranged{curr_ROI}{curr_session} = bsxfun(@minus,df_f_ROI_trial_arranged{curr_ROI}{curr_session},df_f_ROI_trial_arranged{curr_ROI}{curr_session}(:,16));
        end
    end

    % Get trail averaged df/f for each ROI
    for curr_ROI = 1:length(Final_Mask)
        for curr_session = 1:min(11,length(RecICA_Mov))
            df_f_ROI_average{curr_ROI}(curr_session,:) = nanmean(df_f_ROI_trial_arranged{curr_ROI}{curr_session},1);
            % subtract 1st
            df_f_ROI_sub_1st_average{curr_ROI}(curr_session,:) = nanmean(df_f_ROI_sub_1st_trial_arranged{curr_ROI}{curr_session},1);
            % subtract movement onset frame
            df_f_ROI_sub_16th_average{curr_ROI}(curr_session,:) = nanmean(df_f_ROI_sub_16th_trial_arranged{curr_ROI}{curr_session},1);
        end
    end

    % Plot
    colorvalue = colormap(hot);
    figure
    for curr_ROI = 1:length(Final_Mask)
        subplot(4,4,curr_ROI);
        hold on;
        for curr_session = 1:min(11,length(RecICA_Mov))
            plot(df_f_ROI_sub_1st_average{curr_ROI}(curr_session,:),'color',colorvalue(curr_session*5,:));
        end
    end
    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500']);
    if OP
        saveas(gcf,[Initial '_' Animal '_OPExclude.fig']);    
    else
        saveas(gcf,[Initial '_' Animal '.fig']);
    end
    
    colorvalue = colormap(hot);
    figure
    for curr_ROI = 1:length(Final_Mask)
        subplot(4,4,curr_ROI);
        hold on;
        for curr_session = 1:min(11,length(RecICA_Mov))
            plot(df_f_ROI_average{curr_ROI}(curr_session,:),'color',colorvalue(curr_session*5,:));
        end
    end
    if OP
        savefig(gcf,[Initial '_' Animal '_noSub_OPExclude.fig']);
    else
        savefig(gcf,[Initial '_' Animal '_noSub.fig']);
    end
    close(gcf);
    
    close all;
    clear RecICA_Mov;
    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500']);
    if OP
        save([Initial '_' Animal '_ROI_df_f_OPExclude'],'-v7.3');
    else
        save([Initial '_' Animal '_ROI_df_f'],'-v7.3');
    end
end


