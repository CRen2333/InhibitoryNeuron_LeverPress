function [CROPPED]=CM_CROP_MOVIE(STK)
figure;
image_to_Show=STK(:,:,1);
imagesc(image_to_Show);
colormap ('gray')
axis image
h = imrect(gca, [0 20 20 20]);
api = iptgetapi(h);
api.addNewPositionCallback(@(p) title(mat2str(p)));
fcn = makeConstrainToRectFcn('imrect',...
    get(gca,'XLim'),get(gca,'YLim'));
api.setPositionConstraintFcn(fcn);
wait(h)
RECTCOORD=(api.getPosition());
figure,imagesc(image_to_Show);
colormap ('gray')
hold on
RECTCOORD=round (RECTCOORD);
X1=RECTCOORD(1); 
Y1=RECTCOORD(2);
X2=X1+RECTCOORD(3);
Y2=Y1+RECTCOORD(4);
YVECTOR=[Y1,Y1,Y2,Y2,Y1];
XVECTOR=[X1,X2,X2,X1,X1];
ROI_coord=[X1,X2,Y1,Y2];

CROPPED=STK(min(ROI_coord(3:4)):max(ROI_coord(3:4)),min(ROI_coord(1:2)): max(ROI_coord(1:2)),:);

plot (XVECTOR,YVECTOR,'r')
    clear RECTCOORD
    hold off

end