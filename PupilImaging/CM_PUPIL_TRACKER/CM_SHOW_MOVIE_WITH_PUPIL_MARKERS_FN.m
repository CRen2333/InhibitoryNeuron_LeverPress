% CM  20170201

% example 
% CM_SHOW_MOVIE_WITH_PUPIL_MARKERS_FN(CROPPED,hor_mid,ver_mid,hor_pt1,hor_pt2,ver_pt1,ver_pt2)
% CM_SHOW_MOVIE_WITH_PUPIL_MARKERS_FN(CROPPED,pup_hor_mo_18_midline_vector,pupx_ver_mo_18_midline_vector,pup_hor_mo_18_point1_vector,pup_hor_mo_18_point2_vector,pupx_ver_mo_18_point1_vector,pupx_ver_mo_18_point2_vector)

%,,pup_hor_mo_18_point1_vector,pup_hor_mo_18_point2_vector,pupx_ver_mo_18_point1_vector,pupx_ver_mo_18_point2_vector

function CM_SHOW_MOVIE_WITH_PUPIL_MARKERS_FN(CROPPED,hor_mid,ver_mid,hor_pt1,hor_pt2,ver_pt1,ver_pt2)

clear Moviie
figure;
max_x=size(CROPPED,1);
max_y=size(CROPPED,2);
max_z=size(CROPPED,3)-2;
for frame=1:1:max_z;
 %   subplot (1,2,1)
imagesc(CROPPED (:,:,frame)),colormap('gray');axis image;hold on;
caxis([0,220]);
plot(hor_mid(frame),ver_mid(frame),'r*');hold on;
% plot(pup_hor_mo_18_point1_vector(frame),pupx_ver_mo_18_point1_vector(frame),'b*');hold on;
% plot(pup_hor_mo_18_point2_vector(frame),pupx_ver_mo_18_point2_vector(frame),'y*');hold on;

y0=ver_mid(frame);
x0=hor_mid(frame);

x1=hor_pt1(frame);
x2=hor_pt2(frame);
y1=ver_pt1(frame);
y2=ver_pt2(frame);


plot([x1 x2],[y0 y0],'b*');hold on;
plot([x0 x0],[y1 y2],'y*');hold off;
pause (5/100)
%    Moviie(frame)=getframe(gcf); % leaving gcf out crops the frame in the movie.

figure(gcf)
end

% movie2avi(Moviie,'Moviie.avi'); 