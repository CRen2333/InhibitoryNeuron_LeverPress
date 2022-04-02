%% CM 20170215 pupil tracker
% associated code 
% tiffread
% thresh_tool
% CM_CAL_DIAM_pup_MIX
% CM_SHOW_MOVIE_WITH_PUPIL_MARKERS_FN(CROPPED,pup_hor_midline_vector,pup_ver_midline_vector,pup_hor_point1_vector,pup_hor_point2_vector,pup_ver_point1_vector,pup_ver_point2_vector);
%%
[FileName,PathName] = uigetfile('*.avi','Select the movie file');
% v=VideoReader([PathName FileName]);
% while hasFrame(v)
%       images_big = readFrame(v);
% end
moviee=mmread([PathName FileName]);
for i=1:1:moviee.nrFramesTotal
images_big(:,:,i)=moviee.frames(i).cdata(:,:,1);
end
big_time_stamps=1:1:size(images_big,3);
time_per_frame=1/30;
big_time_stamps=(big_time_stamps*time_per_frame)-1/time_per_frame;
rotation_angle=-20; % to adjust
images_big_rot=imrotate(images_big,rotation_angle);figure;imagesc(images_big_rot(:,:,1));axis image;colormap('gray')
%%
[CROPPED]=CM_CROP_MOVIE(images_big_rot); % limit to the eye
[Threshold] = thresh_tool(CROPPED(:,:,1),'gray');
oops=CROPPED<Threshold;

%%
horiz_cross=squeeze(sum(double(oops),1));figure;imagesc((horiz_cross));colormap('jet');
tt=[1 1;size(horiz_cross,2)/8 1;2*(size(horiz_cross,2)/8) 1;3*(size(horiz_cross,2)/8) 1;4*(size(horiz_cross,2)/8) 1;5*(size(horiz_cross,2)/8) 1;6*(size(horiz_cross,2)/8) 1;7*(size(horiz_cross,2)/8) 1;size(horiz_cross,2) 1;size(horiz_cross,2) size(horiz_cross,1)/2; size(horiz_cross,2)  size(horiz_cross,1);size(horiz_cross,2)/2  size(horiz_cross,1);1 size(horiz_cross,1);1 size(horiz_cross,1)/2];
hpoly_hor = impoly(gca,tt);
wait(hpoly_hor);
BWdia_hor = createMask(hpoly_hor);
BWdia2_hor=BWdia_hor;
BWdia2_hor=double(BWdia2_hor);
BWdia2_hor(BWdia_hor~=1)=nan;
image_to_diam_hor=horiz_cross.*BWdia2_hor;
figure; imagesc(image_to_diam_hor);


vert_cross=squeeze(sum(double(oops),2));figure;imagesc(vert_cross);colormap('gray');
tt_vert=[1 1;size(vert_cross,2)/8 1;2*(size(vert_cross,2)/8) 1;3*(size(vert_cross,2)/8) 1;4*(size(vert_cross,2)/8) 1;5*(size(vert_cross,2)/8) 1;6*(size(horiz_cross,2)/8) 1;7*(size(vert_cross,2)/8) 1;size(vert_cross,2) 1;size(vert_cross,2) size(vert_cross,1)/2; size(vert_cross,2)  size(vert_cross,1);size(vert_cross,2)/2  size(vert_cross,1);1 size(vert_cross,1);1 size(vert_cross,1)/2];

hpoly_ver = impoly(gca,tt_vert);
wait(hpoly_ver);
BWdia_ver = createMask(hpoly_ver);
BWdia2_ver=BWdia_ver;
BWdia2_ver=double(BWdia2_ver);
BWdia2_ver(BWdia_ver~=1)=nan;
image_to_diam_ver=vert_cross.*BWdia2_ver;
figure; imagesc(image_to_diam_ver);

FIRST_INDEX=1;
LAST_INDEX=size(image_to_diam_hor,1);
LAST_INDEX_ver=size(image_to_diam_ver,1);

smoothing=1;
threshold_ratio=3;
assignName_hor='pup_hor';
assignName_ver='pup_ver';
time_per_frame=1;
windowSize=2;
windowStep=1;
 %,umPerPix)
 

CM_CAL_DIAM_pup_MIX(FIRST_INDEX,LAST_INDEX,smoothing,threshold_ratio,assignName_hor,time_per_frame,windowSize,windowStep,image_to_diam_hor',1,1,'Horiz') %,umPerPix)
CM_CAL_DIAM_pup_MIX(FIRST_INDEX,LAST_INDEX_ver,smoothing,threshold_ratio,assignName_ver,time_per_frame,windowSize,windowStep,image_to_diam_ver',1,1,'Vert') %,umPerPix)

CM_SHOW_MOVIE_WITH_PUPIL_MARKERS_FN(CROPPED,pup_hor_midline_vector,pup_ver_midline_vector,pup_hor_point1_vector,pup_hor_point2_vector,pup_ver_point1_vector,pup_ver_point2_vector);

%% LID DIAMETER
figure;imagesc(images_big_rot(:,:,1));axis image;colormap('gray')
[EYE_OP]=CM_CROP_MOVIE(images_big_rot);
[Threshold_lid] = thresh_tool(images_big_rot(:,:,1),'gray');

Lid_mov=EYE_OP<Threshold_lid;axis image
figure;subplot(1,2,1);imagesc(Lid_mov(:,:,1));subplot(1,2,2);imagesc(EYE_OP(:,:,1))
vert_cross_LID=squeeze(sum(double(Lid_mov),2));figure;imagesc(vert_cross_LID);colormap('gray');
tt_LID=[1 1;size(vert_cross_LID,2)/8 1;2*(size(vert_cross_LID,2)/8) 1;3*(size(vert_cross_LID,2)/8) 1;4*(size(vert_cross_LID,2)/8) 1;5*(size(vert_cross_LID,2)/8) 1;6*(size(vert_cross_LID,2)/8) 1;7*(size(vert_cross_LID,2)/8) 1;size(vert_cross_LID,2) 1;size(vert_cross_LID,2) size(vert_cross_LID,1)/2; size(vert_cross_LID,2)  size(vert_cross_LID,1);size(vert_cross_LID,2)/2  size(vert_cross_LID,1);1 size(vert_cross_LID,1);1 size(vert_cross_LID,1)/2];
hpoly_ver_LID = impoly(gca,tt_LID);


wait(hpoly_ver_LID);
BWdia_ver_LID = createMask(hpoly_ver_LID);
BWdia2_ver_LID=BWdia_ver_LID;
BWdia2_ver_LID=double(BWdia2_ver_LID);
BWdia2_ver_LID(BWdia_ver_LID~=1)=nan;
image_to_diam_ver_LID=vert_cross_LID.*BWdia2_ver_LID;
figure; imagesc(image_to_diam_ver_LID);

FIRST_INDEX=1;
LAST_INDEX_ver_LID=size(image_to_diam_ver_LID,1);
smoothing=1;
threshold_ratio=8;
assignName_ver_LID='LID_ver';
time_per_frame=1;
windowSize=1;
windowStep=1;
 %,umPerPix)
CM_CAL_DIAM_pup(FIRST_INDEX,LAST_INDEX_ver_LID,smoothing,threshold_ratio,assignName_ver_LID,time_per_frame,windowSize,windowStep,image_to_diam_ver_LID',1,1,'LID') %,umPerPix)




