function []=CM_CAL_DIAM_pup(FIRST_INDEX,LAST_INDEX,smoothing,threshold_ratio,assignName,time_per_frame,windowSize,windowStep,longImageCh1,um_per_pix,DON,string_for_title) %,umPerPix)
nLines = size (longImageCh1,1);    % total number of lines in data;
nLinesPerBlock = round(windowSize / time_per_frame);   % how many lines in each block?
disp (['SM ' num2str(smoothing) '; Thres  ' num2str(threshold_ratio) '; LinesPerBlock ' num2str(nLinesPerBlock) '; time_per_frame ' num2str(time_per_frame*1000) ' ms; timePerBlock ' num2str(nLinesPerBlock*time_per_frame*1000) ' ms'])
windowStartPoints = round(1:windowStep/time_per_frame : nLines-nLinesPerBlock) ;  % where do the windows start (in points?)
nchar=fprintf('ca commence');

for i = 1:length(windowStartPoints)
    if ~mod(i,round(length(windowStartPoints)/50))
        fprintf(repmat('\b',1,nchar))
        string_to_display=['' num2str(round(100*i/length(windowStartPoints))),' percent ', num2str(windowStartPoints(i)) ' lines out of ' num2str(windowStartPoints(end))];
        nchar=fprintf(string_to_display);
    end
    
    w = windowStartPoints(i);         % which line to start this window?
    global_line_start=w;
    global_line_end=w-1+nLinesPerBlock;
    blockData_to_use=longImageCh1 (global_line_start:global_line_end,:);
    blockDataCut = blockData_to_use(:,FIRST_INDEX:LAST_INDEX);
    blockDataMean = mean(blockData_to_use,1);   % take mean of several lines imageCh == 1
    if (i==1)
         CM_big_matrix = zeros (length(windowStartPoints),length(blockDataCut));
        analysisData = 0*windowStartPoints;      % create space to hold data of the variable to return (diameter or velocity)
        point1_vector = 0*windowStartPoints;% create space to hold data for the left limit of the vessel
        point2_vector = 0*windowStartPoints;% create space to hold data for the right limit of the vessel
        intensityData_CM = 0*windowStartPoints; % create space to hold data CM_20121126
    end
    
    blockDataMean = blockDataMean(FIRST_INDEX:LAST_INDEX);  % cut out only portion for this object
    [analysisData(i) point1_vector(i) point2_vector(i) smoothed_profile] = calcFWHM(blockDataMean,smoothing,threshold_ratio);
    intensityData_CM(i) = mean(blockDataMean);
     CM_big_matrix (i,:)=blockDataMean;
    if (i==1)
    CM_big_matrix_smoothed = zeros (length(windowStartPoints),length(smoothed_profile));
    end
    CM_big_matrix_smoothed (i,:)=smoothed_profile;

end

CM_big_matrix_smoothed(1:smoothing,:)=0;
CM_big_matrix_smoothed(end-smoothing:end,:)=0;

%analysisData = analysisData%* umPerPix;     % convert units (currently in pixels) to millivolts

% assignin('base',[assignName '_' 'ch' num2str(imageCh) '_umPerPix'],umPerPix);   %
assignin('base',[assignName '_diameter_pix'],analysisData);   % 
assignin('base',[assignName '_point1_vector'],point1_vector);   % index/sec
assignin('base',[assignName '_point2_vector'],point2_vector);   % index/sec
midline=(point2_vector+point1_vector)/2;
assignin('base',[assignName '_midline_vector'],midline);
assignin('base',[assignName '_Mean_int_mv'],intensityData_CM);
assignin('base',[assignName '_big_matrix'],CM_big_matrix);
assignin('base',[assignName '_big_matrix_SM'],CM_big_matrix_smoothed);

if ~(um_per_pix==1)
    assignin('base',[assignName '_diameter_um'],(analysisData*um_per_pix));   % mv / second
end

% make a time axis that matcheds the diameter info
time_axis = windowSize/2 + windowStep*(0:length(analysisData)-1);
assignin('base',[assignName '_time_axis'],time_axis);

if (DON)
    CM_SHOW_RAW_DATA(CM_big_matrix_smoothed,point1_vector,point2_vector,analysisData,time_axis,string_for_title)
end

disp ' ... DONE'

end


function [width point1 point2 data] = calcFWHM(data,smoothing,threshold_ratio)
% function which takes data and calculates the full-width, half max value
% half-max values are found looking in from the sides, i.e., the program will work
% even if the data dips to a lower value in the middle

point1 = [];
point2 = [];

data = double(data);

% smooth data, if appropriate
if (smoothing==0)
    smoothing = 1;    % smoothing(none)
end
if (smoothing > 1)
    data = conv(data,rectwin(smoothing) ./ smoothing);
end

baseline_to_sub=min(data(smoothing:(length(data)-smoothing)));
data = data-baseline_to_sub; % The minimum is calculated on the second to the avant dernier point
threshold = max(data)/threshold_ratio; % changed from 2 to 4 BY
assignin ('base','ratio_calc_FWHM',threshold_ratio)
%%
aboveI = find(data > threshold);    % all the indices where the data is above half max

if isempty(aboveI)
    % nothing was above threshold!
    width =0;
    point1=0;
    point2=0;
    return
end

firstI = aboveI(1);                 % index of the first point above threshold
lastI = aboveI(end);                % index of the last point above threshold

if (firstI-1 < 1) || (lastI+1) > length(data)
    % interpolation would result in error, set width to zero and just return ...
    width = 0;
    point1=0;
    point2=0;
    return
end

% use linear intepolation to get a more accurate picture of where the max was
% find value difference between the point and the threshold value,
% and scale this by the difference between integer points ...
point1offset = (threshold-data(firstI-1)) / (data(firstI)-data(firstI-1));
point2offset = (threshold-data(lastI)) / (data(lastI+1)-data(lastI));

point1 = firstI-1 + point1offset;
point2 = lastI + point2offset;

width = point2-point1;
%width = lastI-firstI; % CM 20th may

end


function CM_SHOW_RAW_DATA(big_matrix_name,point_1_name,point_2_name,diameter_name,time_name,string_for_title)

figure;
set(gcf,'Name', ['IM_' string_for_title],'Pos',[100 200 700 400],'color','w','NumberTitle','off')
vbig_matrix_name=big_matrix_name;
vpoint_1_name=point_1_name;
vpoint_2_name=point_2_name;
vdiameter_name=diameter_name;
vtime_name=time_name;

subplot(3,1,1)
imagesc (vtime_name,[1:1:size(vbig_matrix_name,2)],vbig_matrix_name');
hold on;
plot(vtime_name,vpoint_1_name,'k');
hold on
plot(vtime_name,vpoint_2_name,'k');

% figure;imagesc (line_2_time_axis,[1:1:size(line_2_ch1_big_matrix_SM,2)],line_2_ch1_big_matrix_SM')

hold off
title (strrep (['IM_' string_for_title],'_','-'));

subplot(3,1,2)
plot(vtime_name,vdiameter_name)
hold on
midline=(vpoint_2_name-vpoint_1_name)/2+vpoint_1_name;
xlabel('Time')
ylabel('Diameter (pix)')
subplot(3,1,3)
plot(vtime_name,(midline-nanmean(midline,2)),'r');
%plot(vtime_name,(midline),'r');

xlabel('Time')
ylabel('Midline (pix)')
ylim ([-60 60])
linkaxesInFigure('x')
end