% Align craniotomy to landmarks on skull to get bregma and coordinates

% Clear
clear all
close all
clc

Initial = 'CR';
Animal = '3526549-O';

clear tformSimilarity movingRegisteredSmilarity
disp(Animal);

cd('Z:\People\Chi\WFLP_IN\Picture_with_Bregma'); 
fixed = imread([Initial '_' Animal '_Bone.jpg']);
fixed = fixed(:,:,2);
moving = imread([Initial '_' Animal '_Clear.jpg']);
moving = moving(:,:,2);
% rotate if necessary
imshow(fixed); % check if roration is needed
% Sometimes the picture needs to be rotated by 90 degrees
% fixed = imrotate(fixed,-90);
% moving = imrotate(moving,-90);
% imshow(fixed);

% rotate fixed to make the midline vertical
fixed_2 = imrotate(fixed,4);
imshow(fixed_2);
imshow(moving);
moving_2 = imrotate(moving,0);
imshow(moving_2);

fixed = fixed_2;
moving = moving_2;
% manually set control points
cpselect(moving,fixed);

tform = fitgeotrans(movingPoints,fixedPoints,'Projective');
% alternative methods
% tform = fitgeotrans(movingPoints,fixedPoints,'Similarity'); 
% tform = fitgeotrans(movingPoints,fixedPoints,'Polynomial',2);
moveing_registered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
figure; % show differences
imshowpair(fixed,moveing_registered,'diff');

[~,ROI_Mask]=CM_CROP_MOVIE(moveing_registered); % get rid of surroundings
[row,col] = find(ROI_Mask);
fixed_cropped = fixed(min(row):max(row),min(col):max(col));
moveing_registered_cropped = moveing_registered(min(row):max(row),min(col):max(col));
imwrite(fixed_cropped,[Initial '_' Animal '_Bone_aligned.jpg']);
imwrite(moveing_registered_cropped,[Initial '_' Animal '_Clear_aligned.jpg']);
% use cropped jpg pictures in Adode Illustrator and overlay coordinate mesh
% on it (_Clear_Marked).
